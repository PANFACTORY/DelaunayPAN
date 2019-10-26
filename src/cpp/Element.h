//*****************************************************************************
//Title		:src/cpp/Element.h
//Author	:Tanabe Yuta
//Date		:2019/10/26
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include "Node.h"


namespace DelaunayPAN{
	class ElementClass{
public:
		ElementClass();
		~ElementClass();

		int node[3];											//���_�̐ߓ_�ԍ�
		int neighbor[3];										//�אڗv�f�ԍ�
		double angle[3];										//���p
		bool side[3];											//true:�ӂ����E		false:�ӂ��񋫊E
		bool active;											//���E�O�v�f������	true�F���E��	false�F���E�O
		bool check;												//���E���O����ς݂��L������	true�F����ς�	false�F������

		void setnode(int _node0, int _node1, int _node2);		//���_�̎�荞��
		void setneighbor(int _neighbor0, int _neighbor1, int _neighbor2);		//�אڗv�f�ԍ��̎�荞��
		void setside(bool _side0, bool _side1, bool _side2);	//���E�ӂ���荞��
		void copy(ElementClass _originalelement);				//�I�u�W�F�N�g���R�s�[����
		void getangle(std::vector<NodeClass> _node);					//���p���v�Z

		int inouton(int _nodenum, std::vector<NodeClass> _node);		//�v�f�Ɠ_�̈ʒu�֌W��Ԃ�
		int oppositenode(int _elementname);						//�w�肳�ꂽ�v�f�Ƃ̈ʒu�֌W��Ԃ�

		double space(std::vector<NodeClass> _node);					//�v�f�ʐς�Ԃ�
		int nodeorder(int _nodenum);							//�w�肳�ꂽ���_�����Ԗڂ��Ԃ�
	};
}