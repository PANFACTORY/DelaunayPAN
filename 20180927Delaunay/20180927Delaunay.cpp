//*****************************************************************************
//Title		:Delaunay
//Purpose	:Create mesh with Delaunay Triangulation Method
//Author	:Tanabe Yuta
//Date		:2018/09/27
//Copyright	:(C) 2018 Tanabe Yuta
//*****************************************************************************
//外部境界：時計回り
//内部境界：反時計回り
//*****************************************************************************

#include "pch.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

#include "ElementClass.h"
#include "NodeClass.h"
#include "BoundaryClass.h"
#include "DelaunayClass.h"

using namespace std;


//*****************************************************************************
//DXFファイル出力
//*****************************************************************************
void exportdxf(vector<ElementClass> _element, vector<NodeClass> _node) {
	double magx = 50.0, magy = 50.0, offsetx = 100.0, offsety = 100.0;
	FILE *fp;
	fopen_s(&fp, "Mesh.dxf", "w");
	fprintf(fp, "0\nSECTION\n2\nHEADER\n9\n$ACADVER\n1\nAC1009\n0\nENDSEC\n0\nSECTION\n2\nTABLES\n0\nENDSEC\n0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n");
	for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int color = 1;
			if (_element[i].side[(j + 2) % 3] == true) {
				color = 7;
			}
			fprintf(fp, "0\nLINE\n8\n0\n62\n%d\n10\n%lf\n20\n%lf\n11\n%lf\n21\n%lf\n", color, _node[_element[i].node[j]].x * magx + offsetx, _node[_element[i].node[j]].y * magy + offsety, _node[_element[i].node[(j + 1) % 3]].x * magx + offsetx, _node[_element[i].node[(j + 1) % 3]].y * magy + offsety);
		}
	}
	fprintf(fp, "0\nENDSEC\n0\nEOF");
	fclose(fp);
}


//*****************************************************************************
//要素の出力
//*****************************************************************************
void showelements(vector<ElementClass> _element, vector<NodeClass> _node) {
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
	for (int i = 0; i < _element.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int linename = 1;
			if (_element[i].side[(j + 2) % 3] == true) {
				linename = 2;
			}
			fprintf(fp, "set arrow from %lf,%lf to %lf,%lf arrowstyle %d\n", _node[_element[i].node[j]].x, _node[_element[i].node[j]].y, _node[_element[i].node[(j + 1) % 3]].x, _node[_element[i].node[(j + 1) % 3]].y, linename);
		}
	}
	fprintf(fp, "plot '-' with points pointtype 1\n");
	fprintf(fp, "e\n");

	fprintf(fp, "set nomultiplot\n"); // マルチプロットモード終了
	fprintf(fp, "exit\n"); // gnuplotの終了
	fflush(fp);
}


//*****************************************************************************
//節点の取り込み
//*****************************************************************************
void importnode(vector<NodeClass> &_node, string _fname) {
	string tmp;

	ifstream fin(_fname);
	if (!fin) {
		cout << "Node File Open Error\n";
	}

	while (getline(fin, tmp)) {
		istringstream ssi(tmp);
		string tmpx, tmpy;
		NodeClass tmpnode;

		ssi >> tmpx >> tmpy;
		tmpnode.x = stod(tmpx);
		tmpnode.y = stod(tmpy);
		_node.push_back(tmpnode);
	}
	fin.close();
}


//*****************************************************************************
//境界の取り込み
//*****************************************************************************
void importboundary(vector<BoundaryClass> &_boundary, string _fname, bool _type) {
	BoundaryClass tmpboundary;
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
	vector<NodeClass> node;							//節点
	vector<ElementClass> element;					//要素
	vector<BoundaryClass> boundary;					//境界
	DelaunayClass mesher;							//DelaunayTriangulationによるメッシャー

	//----------節点座標の取り込み----------
	importnode(node, "Model2/nodes.dat");

	//----------境界の取り込み----------
	importboundary(boundary, "Model2/externalboundary0.dat", true);
	//importboundary(boundary, "Model2/externalboundary1.dat", false);
	//importboundary(boundary, "Model1/internalboundary1.dat", false);

	//----------SuperTriangleの生成----------
	mesher.delaunaymain(node, element, boundary, 0.5);
	
	//----------結果の表示----------
	cout << element.size() << "\n";
	showelements(element, node);
	//exportdxf(element, node);
}