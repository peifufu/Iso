#include"Jie_Ds.h"
#include"Vec_nD.h"
#include<iostream>
using namespace std;
int main()
{
	class Jie_Ds Jie_test;
	vector<Data> Convexhull, ConvexhullTest,ConvexhullExp;
	//读入井位点信息
	Convexhull = Jie_test.LoadModel("data.txt");
	std::ofstream out1("point_ori.obj");
	for (int j = 0; j < Convexhull.size(); j++)
	{
		Data qq;
		qq = Convexhull[j];
		out1 << "v " << fixed << setprecision(5)<<qq.X << " " << fixed << setprecision(5)<<qq.Y << " " << "0" << endl;
	}
	out1.close();
	//计算凸包(凸包计算完成后，凸包上面的点会沿着顺时针、逆时针两种分布，此处是逆时针)
	ConvexhullTest = Jie_test.Withershins(Convexhull);
	std::ofstream out("convexhull_ori.obj");
	for (int j = 0; j < ConvexhullTest.size(); j++)
	{
		Data qq;
		qq = ConvexhullTest[j];
		out << "v " << fixed << setprecision(5)<<qq.X << " " << fixed << setprecision(5)<<qq.Y << " " << "0" << endl;
	}
	int count_ori = 1;
	for (int j = 0; j < ConvexhullTest.size() - 1; j++)
	{
		Data qq1, qq2;
		qq1 = ConvexhullTest[j];
		qq2 = ConvexhullTest[j + 1];
		out << "l " << count_ori << " " << count_ori + 1 << endl;
		count_ori++;
	}
	out.close();

	double Sum_x=0,Sum_y=0,Aver_x=0, Aver_y=0;
	for (int j = 0; j < ConvexhullTest.size()-1; j++)
	{
		Sum_x += ConvexhullTest[j].X;
		Sum_y += ConvexhullTest[j].Y;
	}
	Aver_x = Sum_x / (ConvexhullTest.size()-1);
	Aver_y = Sum_y / (ConvexhullTest.size()-1);
	std::ofstream out9("center.obj");
	out9 << "v " << fixed << setprecision(5) << Aver_x << " " << fixed << setprecision(5) << Aver_y << " " << "0" << endl;
	out9.close();



	for (int i = 0; i < ConvexhullTest.size(); i++)
	{
		Data cp;
		cp.X = ConvexhullTest[i].X;
		cp.Y = ConvexhullTest[i].Y;
		double d_v_x = ConvexhullTest[i].X - Aver_x;
		double d_v_y = ConvexhullTest[i].Y - Aver_y;
		double d_n_v_x = d_v_x / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		double d_n_v_y = d_v_y / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		cp.X += 3* d_n_v_x;
		cp.Y += 3* d_n_v_y;
		ConvexhullExp.push_back(cp);
	}

	std::ofstream out12("convexhull_exp.obj");
	for (int j = 0; j < ConvexhullExp.size(); j++)
	{
		Data qq;
		qq = ConvexhullExp[j];
		out12 << "v " << fixed << setprecision(5) << qq.X << " " << fixed << setprecision(5) << qq.Y << " " << "0" << endl;
	}
	int count_exp = 1;
	for (int j = 0; j < ConvexhullExp.size() - 1; j++)
	{;
		out12 << "l " << count_exp << " " << count_exp + 1 << endl;
		count_exp++;
	}
	out12.close();



	
	//cout << "hello world!" << endl;
	cout << "before oridata:" << Jie_test.m_oriData.size() << endl;
	cout << "before border:"<<Jie_test.m_Border.size() << endl;
	cout << "before oriborder:" << Jie_test.m_OriBoder.size() << endl;
	cout << "before grid:"<<Jie_test.m_GridPoint.size() << endl;
	

	for (int i = 0; i < Convexhull.size(); i++)
		Jie_test.AddData(Convexhull[i]);

	Jie_test.Charact();
	//Jie_test.CalcBorder();
	//cout << Jie_test.m_Border.size() << endl;
	Jie_test.m_XNum = 200;
	Jie_test.m_YNum = 200;
	Jie_test.SetGridXY();
	std::ofstream out5("convexhull_border.obj");
	for (int i = 0; i < Jie_test.m_Border.size(); i++)
	{
		out5 << "v " << fixed << setprecision(5) << Jie_test.m_Border[i].X << " " << fixed << setprecision(5) << Jie_test.m_Border[i].Y << " " << "0" << endl;
	}
	int count_border = 1;
	for (int i = 0; i < Jie_test.m_Border.size() - 1; i++)
	{
		out5 << "l " << count_border << " " << count_border + 1 << endl;
		count_border += 1;
	}
	out5.close();





	//插值之前网格点的属性值
	/*for (int i = 1; i < Jie_test.m_GridPoint.size() - 2; i++)
	{
	for (int j = 1; j < Jie_test.m_GridPoint.size() - 2; j++)
	{
		cout << Jie_test.m_GridPoint[i][j].Z <<" ";
	}
	cout << endl;
	}*/
	Jie_test.EvaluateNoFault();

	cout << "after oridata:" << Jie_test.m_oriData.size() << endl;
	cout << "after border:" << Jie_test.m_Border.size() << endl;
	cout << "sfter oriborder:" << Jie_test.m_OriBoder.size() << endl;
	cout << "after grid:"<<Jie_test.m_GridPoint.size() << endl;


	//cout << "v " << Jie_test.m_GridPoint[1][1].X << " " << Jie_test.m_GridPoint[1][1].Y << " " << "0" << endl;
	//cout << "v " << Jie_test.m_GridPoint[1][2].X << " " << Jie_test.m_GridPoint[1][2].Y << " " << "0" << endl;
	std::ofstream out2("Grid.obj");
	int count_grid = 1;
	// 按照Grid.h中的计算，横纵行列都是从1开始,最后两行/列也是附加的，都是0,0,0;不要了
	// Grid.h中的m_GridPoint，一级下标对应的是网格的列数据，二级下标对应的是行数据
	for (int i = 1; i < Jie_test.m_GridPoint.size()-2; i++)
	{
		for (int j = 1; j < Jie_test.m_GridPoint.size()-2; j++)
		{
			out2 << "v " << fixed << setprecision(5)<<Jie_test.m_GridPoint[i][j].X << " " << fixed << setprecision(5) << Jie_test.m_GridPoint[i][j].Y << " " << "0" << endl;
			//cout << Jie_test.m_GridPoint[i][j].X << " " << Jie_test.m_GridPoint[i][j].Y << " " << Jie_test.m_GridPoint[i][j].Z << endl;
		}
	}
	int current = Jie_test.m_GridPoint.size() - 3;
	int flag;
	for (int i = 1; i < Jie_test.m_GridPoint.size() - 2; i++)
	{
		flag = 1;
		int j = (i-1)*current+1;
		while (flag < current)
		{
			out2 << "l " << j << " " << j + 1 << endl;
			j += 1;
			flag++;
		}
	}

	
	for (int i = 1; i < Jie_test.m_GridPoint.size() - 2; i++)
	{
		flag = 1;
		int j = i;
		while (flag < current)
		{
			out2 << "l " << j << " " << j + current << endl;
			j += current;
			flag++;
		}
	}
	out2.close();

	//插值之后网格点的属性值
	/*for (int i = 1; i < Jie_test.m_GridPoint.size() - 2; i++)
	{
		for (int j = 1; j < Jie_test.m_GridPoint.size() - 2; j++)
		{
			cout << Jie_test.m_GridPoint[i][j].Z<<" " ;
		}
		cout << endl;
	}*/

	Jie_test.m_Show_MinValue = Jie_test.m_ZMin;//数据点属性值的最小值
	Jie_test.m_Show_MaxValue = Jie_test.m_ZMax;//数据点属性值的最大值
	Jie_test.m_Show_JianGeValue = abs(Jie_test.m_Show_MaxValue - Jie_test.m_Show_MinValue)*0.2;//这里应该是等值线的间隔
	cout << endl;
	cout << endl;
	//cout << Jie_test.m_Show_MaxValue << endl;
	//cout << Jie_test.m_Show_MinValue << endl;
	//cout << Jie_test.m_Show_JianGeValue << endl;
	Jie_test.CalcSameArray();
	//cout << Jie_test.m_Show_MaxValue << endl;
	//cout << Jie_test.m_Show_MinValue << endl;
	//cout << Jie_test.m_Show_JianGeValue << endl;
	vector<double> vd;
	double m_nDenseNum = 1;
	double dvalue = (Jie_test.m_Show_JianGeValue) / (m_nDenseNum + 1);
	double dmax = (Jie_test.m_Show_MaxValue);
	double dmin = (Jie_test.m_Show_MinValue);

	while (dmin< dmax)
	{
		//cout << "dmin: " << dmin << " ";
		vd.push_back(dmin);
		dmin += dvalue;
	}
	//cout << vd.size() << endl;
	Jie_test.SetTrackValue(vd);
	Jie_test.IsolineTracking();
	

	//Jie_test.CreateIsoLine();
	//cout << Jie_test.m_IsoLine.size() << endl;
	//cout << Jie_test.m_lsoLines.size() << endl;

	cout << "======================================================" << endl;
	cout << "以下是关于等值线的一些信息" << endl;
	cout << "等值线的数量：" << Jie_test.m_IsoRealLine.size() << endl;
	cout << "各条等值线的高程值：" << endl;
	for (int i = 0; i < Jie_test.m_IsoRealLine.size(); i++)
	{
		cout << "value = :" << Jie_test.m_IsoRealLine[i].Value << endl;
	}
	//std::cout << Jie_test.Jie_IsoLine.size() << endl;
	std::ofstream out6("isoline_m_IsoLine.obj");
	for (int i = 0; i < Jie_test.m_IsoLine.size(); i++)
	{
		for (int j = 0; j < Jie_test.m_IsoLine[i].Logic.size(); j++)
		{
			out6 << "v " << Jie_test.m_IsoLine[i].Logic[j].X << " " << Jie_test.m_IsoLine[i].Logic[j].Y << " " << "0" << endl;
		}
	}
	out6.close();
	std::ofstream out7("isoline_m_IsoRealLine.obj");
	std::ofstream out10("startp.obj");
	std::ofstream out11("endp.obj");
	int count_realline = 1;
	for (int i = 0; i < Jie_test.m_IsoRealLine.size(); i++)
	{
		//cout << Jie_test.m_IsoRealLine[i].Value << " ";
		out10 << "v " << Jie_test.m_IsoRealLine[i].Logic[0].X << " " << Jie_test.m_IsoRealLine[i].Logic[0].Y << " " << "0" << endl;
		out11 << "v " << Jie_test.m_IsoRealLine[i].Logic[Jie_test.m_IsoRealLine[i].Logic.size()-1].X << " " << Jie_test.m_IsoRealLine[i].Logic[Jie_test.m_IsoRealLine[i].Logic.size()-1].Y << " " << "0" << endl;
		for (int j = 0; j < Jie_test.m_IsoRealLine[i].Logic.size(); j++)
		{
			out7 << "v " << fixed << setprecision(5) << Jie_test.m_IsoRealLine[i].Logic[j].X << " " << fixed << setprecision(5) << Jie_test.m_IsoRealLine[i].Logic[j].Y << " " << "0" << endl;
		}
		for (int j = 0; j < Jie_test.m_IsoRealLine[i].Logic.size() - 1; j++)
		{
			out7 << "l " << count_realline << " " << count_realline + 1 << endl;
			count_realline++;
		}
		count_realline++;
	}
	out7.close();
	out10.close();
	out11.close();

	std::ofstream out3("Equi_point.obj");
	for (int i = 0; i < Jie_test.Jie_IsoLine.size(); i++)
	{
		for (int j = 0; j < Jie_test.Jie_IsoLine[i].Logic.size(); j++)
		{	
			out3 << "v " << fixed << setprecision(5) << Jie_test.Jie_IsoLine[i].Logic[j].X << " " << fixed << setprecision(5) << Jie_test.Jie_IsoLine[i].Logic[j].Y <<" "<<"0"<< endl;
		}
	}
	out3.close();



	//cout << "m_XMin: " << fixed << setprecision(5) << Jie_test.m_XMin << "m_YMin: " << fixed << setprecision(5) << Jie_test.m_YMin << endl;














	cout << endl;
	cout << "=======================================================" << endl;
	cout << "=======================================================" << endl;
	cout << "下面进行的是边界对于等值线的切割操作" << endl;
	vector<TWOVALUE> Jie_border;
	for (int j = 0; j < Convexhull.size(); j++)
	{
		Data qq;
		qq = Convexhull[j];
		TWOVALUE pp;
		pp.X = qq.X;
		pp.Y = qq.Y;
		Jie_border.push_back(pp);
	}
	
	Jie_test.SetOriBoder(Jie_border);
	Jie_test.SetOriISOLine(Jie_test.m_IsoRealLine);

	Jie_test.DleaIso();

	std::ofstream out8("aftercut_isoline.obj");
	for (int i = 0; i < Jie_test.Jie_IsoRealLine.size(); i++)
	{
		for (int j = 0; j < Jie_test.Jie_IsoRealLine[i].Logic.size(); j++)
			out8 << "v " << Jie_test.Jie_IsoRealLine[i].Logic[j].X << " " << Jie_test.Jie_IsoRealLine[i].Logic[j].Y << " " << "0" << endl;
	}
	out8.close();







	getchar();
	return 0;
}