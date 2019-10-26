//*****************************************************************************
//Title		:src/cpp/main.cpp
//Purpose	:Create mesh with Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27～
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
//注記：
//外部境界：時計回りで定義，trueに設定
//内部境界：反時計回りで定義，falseに設定
//
//手順としては
//①節点（Node）と境界（Boundary）を定義
//②要素配列を定義（Element）
//③Mesherオブジェクト（Delaunay）を生成し実行
//	なおDelaunayのdelaunaymainの最後の引数は分割後の要素の辺長さ最大値
//*****************************************************************************
//要修正：
//要素の細かさ・密度分布制御->面積基準？八分木法
//要素の歪み修正->Laplacian法
//*****************************************************************************


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>


#include "Element.h"
#include "Node.h"
#include "Boundary.h"
#include "Delaunay.h"


using namespace std;
using namespace DelaunayPAN;


//*****************************************************************************
//DXFファイル出力
//*****************************************************************************
void exportdxf(vector<Element<double> > _element, vector<Node<double> > _node) {
	double magx = 50.0, magy = 50.0, offsetx = 100.0, offsety = 100.0;
	
    ofstream fout("Mesh.dxf");

    fout << "0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1009\n0\nENDSEC\n0\nSECTION\n2\nTABLES\n0\nENDSEC\n0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n";
    for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int color = 1;
			if (_element[i].side[(j + 2) % 3] == true) {
				color = 7;
			}
			fout << "0\nLINE\n8\n0\n62\n" << color << "\n10\n" << _node[_element[i].node[j]].x * magx + offsetx << "\n20\n" << _node[_element[i].node[j]].y * magy + offsety << "\n11\n" << _node[_element[i].node[(j + 1) % 3]].x * magx + offsetx << "\n21\n" << _node[_element[i].node[(j + 1) % 3]].y * magy + offsety << endl;
		}
	}
    fout << "0\nENDSEC\n0\nEOF";

    fout.close();    
}


//*****************************************************************************
//要素の出力
//*****************************************************************************
void showelements(vector<Element<double> > _element, vector<Node<double> > _node) {
	FILE *fp;
	fp = _popen("gnuplot -persist", "w");							// パイプを開き、gnuplotの立ち上げ
	fprintf(fp, "set title 'Triangle Elements'\n");					//グラフタイトル
	fprintf(fp, "set multiplot\n");									// マルチプロットモード
	fprintf(fp, "set xrange [%lf:%lf]\n", -1.0, 6.0);	     		// 範囲の指定
	fprintf(fp, "set yrange [%lf:%lf]\n", -1.0, 2.0);
	fprintf(fp, "set xlabel \"x\"\n");							    // ラベル表示
	fprintf(fp, "set ylabel \"y\"\n");
	fprintf(fp, "set style arrow 1 nohead linecolor rgb 'red'\n");
	fprintf(fp, "set style arrow 2 nohead linecolor rgb 'black'\n");
	fprintf(fp, "set size ratio - 1\n");
	// 点のプロット
	//境界線描画
	fprintf(fp, "plot '-' w l lw 1 lt rgb 'black'\n");
	for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			if (_element[i].side[(j + 2) % 3] == true) {
				fprintf(fp, "%lf %lf \n %lf %lf \n \n",
					_node[_element[i].node[j]].x, _node[_element[i].node[j]].y, _node[_element[i].node[(j + 1) % 3]].x, _node[_element[i].node[(j + 1) % 3]].y);

			}
		}
	}
	fprintf(fp, "e\n");
	
	//内部要素描画
	fprintf(fp, "plot '-' w l lw 1 lt rgb 'red'\n");
	for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			if (_element[i].side[(j + 2) % 3] == false) {
				fprintf(fp, "%lf %lf \n %lf %lf \n \n",
					_node[_element[i].node[j]].x, _node[_element[i].node[j]].y, _node[_element[i].node[(j + 1) % 3]].x, _node[_element[i].node[(j + 1) % 3]].y);

			}
		}
	}

	fprintf(fp, "e\n");
	fprintf(fp, "set nomultiplot\n"); // マルチプロットモード終了
	fflush(fp);
}


//*****************************************************************************
//節点の取り込み
//*****************************************************************************
void importnode(vector<Node<double> > &_node, string _fname) {
	string tmp;

	ifstream fin(_fname);
	if (!fin) {
		cout << "Node File Open Error\n";
	}

	while (getline(fin, tmp)) {
		istringstream ssi(tmp);
		string tmpx, tmpy;
		ssi >> tmpx >> tmpy;
		_node.push_back(Node<double>(stod(tmpx), stod(tmpy)));
	}
	fin.close();
}


//*****************************************************************************
//境界の取り込み
//*****************************************************************************
void importboundary(vector<Boundary> &_boundary, string _fname, bool _type) {
	Boundary tmpboundary;
	string tmp;

	ifstream fin(_fname);
	if (!fin) {
		cout << _boundary.size() << " Boundary File Open Error\n";
	}

	while (getline(fin, tmp)) {
		tmpboundary.nodelist.push_back(stoi(tmp));
	}
		
	fin.close();
	tmpboundary.type = _type;
	_boundary.push_back(tmpboundary);
}


//*****************************************************************************
//メイン処理
//*****************************************************************************
int main() {
	//----------節点座標の取り込み----------
	vector<Node<double> > node;							//節点
	importnode(node, "sample/Model1/nodes.dat");

	//----------境界の取り込み----------
	vector<Boundary> boundary;					//境界
	importboundary(boundary, "sample/Model1/externalboundary0.dat", true);
	importboundary(boundary, "sample/Model1/externalboundary1.dat", true);
	importboundary(boundary, "sample/Model1/externalboundary2.dat", true);
	importboundary(boundary, "sample/Model1/internalboundary0.dat", false);
	importboundary(boundary, "sample/Model1/internalboundary1.dat", false);

	//----------DelaunayTriangulationの実行----------
	vector<Element<double> > element;					//要素
	clock_t ts = clock();
	Delaunay<double> mesher;							//Mesherオブジェクト
	mesher.delaunaymain(node, element, boundary, 0.05, 1000);
	clock_t te = clock();
	cout << "\ntime cost:\t" << (double)(te - ts) / CLOCKS_PER_SEC << "sec.\n";
	//----------結果の表示----------
	cout << element.size() << "\n";
	showelements(element, node);
	exportdxf(element, node);
}