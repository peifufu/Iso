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
	double x; //横坐标
	double y; //纵坐标
	double z; // 属性值
};
struct HwIsoLine // 等值线线结构
{
	double zValue; // 等值线属性值
	HwIsoPt* pts; // 等值线平面坐标
	int ptn; //坐标个数
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
//等值线数据体
struct IsoLine
{
	TWOVALUE direction0;	//歼灭井位置，从歼灭井到非歼灭井  20140313
	TWOVALUE direction1;	//非歼灭井，从歼灭井到非歼灭井  20140313
	int  Index;				//跟direction0，direction1对应的点在等值线中的位置（用于等值区域符号填充及尖灭方向判断）20130314
	double Value;			//等值线的属性值
	vector<TWOVALUE> Logic; //等值线的逻辑坐标
};
struct Dis //此数据结构代表什么
{
	double dis;//长度
	double per;//百分比
	double angle;//角度
};

class Data
{
public:
	Data(double x = 0, double y = 0, double z = 0) { X = x; Y = y; Z = z; mark = 0; };
	double X;
	double Y;
	double Z;
	int mark;//判断正断层还是逆断层  正断层为1 逆断层为0
	void Init(double x = 0, double y = 0, double z = 0) { X = x; Y = y; Z = z; mark = 0; }

};
class Jie_Ds
{
public:
	vector<Data> m_oriData; //原始数据，井位信息；输入数据
							//vector<FaultPoint> m_FaultData; //断层信息20130814
	vector<vector<THRVALUE>> m_GridPoint;//生成的网格
	int m_XNum; //横向网格数量
	int m_YNum; //纵向网格数量
	int m_DIV; //等值线间隔
	double m_XMin; //x方向最小值
	double m_XMax;  //X方向最大值
	double m_YMin; //y方向最小值
	double m_YMax; //Y方向最大值
	double m_ZMin; //属性最小值
	double m_ZMax; //属性最大值
	double m_Show_MinValue;
	double m_Show_MaxValue;
	double m_Show_JianGeValue;
	vector<Data> m_Border;
	vector<Data> m_OriBoder;					//原始的真实边界（可能是凹断层）20140801
	bool m_IsConvex;							//原始的边界是否为凸包20140801
	bool m_IsK;								//是否使用克里金估值 20150323

	bool m_bUseFault;
	double m_valuedis;//等值线间隔 20140724
	double m_B;		 //直线型变异函数B值20140728
	vector<Data> m_suV;
	vector<double> m_ni;
	vector<IsoLine>m_IsoLine;
	//存储所有的等值线,每条等值线存储实际值，主要是等值点在断层上不好使用虚拟存储
	vector<IsoLine>m_IsoRealLine;
	vector<double> m_TrackValue;				//要寻找的等值线值
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
