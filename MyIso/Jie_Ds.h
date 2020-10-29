#pragma once
#include "afxtempl.h"
#include<iomanip>
#include"Vec_nD.h"
#include<vector>
#include<string>
#include<fstream>
#include <sstream>
#include<cstdio>
#include<ctime>
#include<iostream>
using namespace std;
#define random(x) rand()%(x)
const int maxn = 1e4 + 20;
struct point
{
	float x, y;
	point() {}
	point(float xx, float yy)
		:x(xx), y(yy) {}
};

struct HwIsoPt
{
	double x; //������
	double y; //������
	double z; // ����ֵ
};
struct HwIsoLine // ��ֵ���߽ṹ
{
	double zValue; // ��ֵ������ֵ
	HwIsoPt* pts; // ��ֵ��ƽ������
	int ptn; //�������
};
struct TWOVALUE
{
	double X;
	double Y;
};
struct THRVALUE
{
	double X;
	double Y;
	double Z;
};
//��ֵ��������
struct IsoLine
{
	TWOVALUE direction0;	//����λ�ã��Ӽ��𾮵��Ǽ���  20140313
	TWOVALUE direction1;	//�Ǽ��𾮣��Ӽ��𾮵��Ǽ���  20140313
	int  Index;				//��direction0��direction1��Ӧ�ĵ��ڵ�ֵ���е�λ�ã����ڵ�ֵ���������估�������жϣ�20130314
	double Value;			//��ֵ�ߵ�����ֵ
	vector<TWOVALUE> Logic; //��ֵ�ߵ��߼�����
};
struct Dis //�����ݽṹ����ʲô
{
	double dis;//����
	double per;//�ٷֱ�
	double angle;//�Ƕ�
};

class Data
{
public:
	Data(double x = 0, double y = 0, double z = 0) { X = x; Y = y; Z = z; mark = 0; };
	double X;
	double Y;
	double Z;
	int mark;//�ж����ϲ㻹����ϲ�  ���ϲ�Ϊ1 ��ϲ�Ϊ0
	void Init(double x = 0, double y = 0, double z = 0) { X = x; Y = y; Z = z; mark = 0; }

};
class Jie_Ds
{
public:
	vector<Data> m_oriData; //ԭʼ���ݣ���λ��Ϣ����������
							//vector<FaultPoint> m_FaultData; //�ϲ���Ϣ20130814
	vector<vector<THRVALUE>> m_GridPoint;//���ɵ�����
	int m_XNum; //������������
	int m_YNum; //������������
	int m_DIV; //��ֵ�߼��
	double m_XMin; //x������Сֵ
	double m_XMax;  //X�������ֵ
	double m_YMin; //y������Сֵ
	double m_YMax; //Y�������ֵ
	double m_ZMin; //������Сֵ
	double m_ZMax; //�������ֵ
	double m_Show_MinValue;
	double m_Show_MaxValue;
	double m_Show_JianGeValue;
	vector<Data> m_Border;
	vector<Data> m_OriBoder;					//ԭʼ����ʵ�߽磨�����ǰ��ϲ㣩20140801
	bool m_IsConvex;							//ԭʼ�ı߽��Ƿ�Ϊ͹��20140801
	bool m_IsK;								//�Ƿ�ʹ�ÿ�����ֵ 20150323

	bool m_bUseFault;
	double m_valuedis;//��ֵ�߼�� 20140724
	double m_B;		 //ֱ���ͱ��캯��Bֵ20140728
	vector<Data> m_suV;
	vector<double> m_ni;
	vector<IsoLine>m_IsoLine;
	//�洢���еĵ�ֵ��,ÿ����ֵ�ߴ洢ʵ��ֵ����Ҫ�ǵ�ֵ���ڶϲ��ϲ���ʹ������洢
	vector<IsoLine>m_IsoRealLine;
	vector<double> m_TrackValue;				//ҪѰ�ҵĵ�ֵ��ֵ
	CArray<double,double> m_sameArray;
	vector<HwIsoLine> m_lsoLines;
	vector<IsoLine> Jie_IsoLine;
	vector<TWOVALUE>Jie_OriBoder;
	vector<IsoLine>Jie_IsoRealLine;
	vector<double> Jie_RNum;
	vector<Vector2D> originPoint;
	vector<Vector2D> curvePoint;
	vector<point> input_vertice;
	vector<point> control_point;
	vector<point> b_spline;
public:
	void GetMinD(vector<double> CosS, vector<double> DisT, int &index);
	vector<Data> Withershins(vector<Data> m_point);
	vector<Data> LoadModel(const char* sFileName);
	void Charact();
	void DataOpt();
	void DataRec();
	void CalcBorder();
	void AddData(Data &t);
	double GetDis(double x1, double y1, double x2, double y2);
	void AddPt(vector<Data> &convexBag);
	double Angle(Data &p0, const Data &p1, const Data &p2);
	double bezier3funcX(double uu, Vector2D *controlP);
	double bezier3funcY(double uu, Vector2D *controlP);
	void createCurve();
	void deCasteljau();
	void deBoor();
	void OptimizeBoder(vector<Data> &convexBag,double e);
	void OptimizeBorderBezier(vector<Data> &convexBag,double e);
	void DividedFourParts(int n, vector<int>& FourPort);
	void BordersPointDis(vector<int>& FourPart, vector<Data>& convexBag, vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
		vector<Dis>& upDis, vector<Dis>&downDis, vector<Dis>&leftDis, vector<Dis>&rightDis);
	void BordersChar(vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
		vector<Dis>& upDis, vector<Dis>&downdis, vector<Dis>&leftDis, vector<Dis>&rightDis);
	void BordersPoints(vector<int>& FourPart, vector<Data>& convexBag, vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
		vector<Dis>& upDis, vector<Dis>&downdis, vector<Dis>&leftDis, vector<Dis>&rightDis,
		vector<Data>& DL, vector<Data>& DR, vector<Data>&DD, vector<Data>& DU,
		vector<TWOVALUE>& sim1, vector<TWOVALUE>& sim2, vector<TWOVALUE>& sim3, vector<TWOVALUE>& sim4, int XM, int YM);
	bool Inv(vector<vector<double>>&M);
	double InsertEst(vector<Data>& suV, TWOVALUE& D, vector<double>& ni);
	void PreMatrix(vector<Data>& suV, vector<double>& ni);
	double DisInv(TWOVALUE D);
	void EvaluateNoFault();
	void SetGridXY();
	double GetMagnitude(double fNumber);
	float FindStep(float StepMin, bool bUporDown);
	void CalcSameArray();
	void SetTrackValue(vector<double> Track);
	void TrackRight(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	void TrackLeft(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	void TrackDown(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	void TrackUp(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	void GetReallyPoint(TWOVALUE A, TWOVALUE &B);
	void TrackY(TWOVALUE A, double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	void TrackX(TWOVALUE A, double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used);
	//void TrackPoint(double Value, IsoLine &Line, IsoLine &IsoReal, double flag_x[][1000], double flag_y[][1000], TWOVALUE First);
	void TrackPoint(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used, TWOVALUE First);
	void SignBorder(vector<vector<double>>& biaoji1, vector<vector<double>>& biaoji2);
	void EquivalentPoints(double Value, vector<TWOVALUE>&VirtualIJ);
	void EquivalentPoints(double Value, vector<THRVALUE>&VirtualIJ);
	void TrackOneValue(double Value);
	void IsolineTracking();
	void CreateIsoLine();
	void SetOriBoder(vector<TWOVALUE> OriBoder);
	void SetOriISOLine(vector<IsoLine> IsoRealLine);
	bool L2L_Intersect(TWOVALUE p1_s, TWOVALUE p1_e, TWOVALUE p2_s, TWOVALUE p2_e, TWOVALUE &Point);
	bool GetCrossPt(TWOVALUE Star, TWOVALUE End, vector<TWOVALUE> Line, TWOVALUE &A);
	bool ISIntersect(TWOVALUE p1_s, TWOVALUE p1_e, TWOVALUE p2_s, TWOVALUE p2_e, double &index);
	bool IsinLineK(TWOVALUE Star, TWOVALUE End, TWOVALUE Point);
	bool IsInside(TWOVALUE A, vector<TWOVALUE> Line);
	void DleaIso(IsoLine &OneIso, vector<IsoLine> &NewIso);
	void DleaIso();
	vector<double> Randnum(int k);
};
