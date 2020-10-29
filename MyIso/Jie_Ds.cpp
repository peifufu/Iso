#include "Jie_Ds.h"

void Jie_Ds::GetMinD(vector<double> CosS, vector<double> DisT, int &index)
{
	double MinCos = CosS[0];
	index = 0;
	for (int i = 0; i < (int)CosS.size(); i++)
	{
		if (CosS[i] < MinCos)
		{
			MinCos = CosS[i];
			index = i;
		}
		else if (CosS[i] == MinCos && DisT[i] >= DisT[index])
		{
			index = i;
		}
	}
}
vector<Data> Jie_Ds::Withershins(vector<Data> m_point)
{

	Data Temppoint;
	int Numindex = (int)m_point.size();
	int Num = (int)m_point.size();
	Temppoint = m_point[0];
	vector<Data> Point;
	int StarIndex = 0;
	//选择最小的x,最大的y作为初始点,存储在m_point[0]中
	for (int i = 1; i< (int)m_point.size(); i++)
	{
		if (Temppoint.X == m_point[i].X && Temppoint.Y < m_point[i].Y)
		{
			Temppoint = m_point[i];
			StarIndex = i;
		}
		else if (Temppoint.X > m_point[i].X)
		{
			Temppoint = m_point[i];
			StarIndex = i;
		}
	}
	//计算所有点和该点组成的向量和初始向量（pt = ptStartA-ptEndA，且使得该向量的x值为正值）之间的的夹角，取一个最大的夹角(最小的cos x = (a*b)/(|a|*|b|))
	//这里的初始向量指的应该是取一个y轴正方向的向量即可，比如说（0,1）
	Point.push_back(Temppoint);
	m_point[StarIndex].mark = 1;	//被使用过
	if ((int)m_point.size() == 1)
	{
		return Point;		//处理单个井位估值 20131108
	}
	vector<double> CosS, NiCosS, Dis, NiDis;
	vector<int> IndexS, NiIndexS;
	//找到第二个点，将其放到m_point[1];
	for (int i = 0; i< Num; i++)
	{
		if (m_point[i].mark == 1)
		{
			continue;	//被使用过
		}
		double x1 = 0.0;
		double y1 = 1.0;
		double x2 = m_point[i].X - Point[0].X;
		double y2 = m_point[i].Y - Point[0].Y;
		double d = -x2;		//叉乘结果，这里主要根据符号来判断(x2,y2)所构成的向量在(x1,y1)所构成向量的顺时针还是逆时针方向
							//若小于0，则在顺时针方向，若大于0，则在逆时针方向，若等于0，则共线

		double d1 = sqrt(pow(x1, 2) + pow(y1, 2));
		double d2 = sqrt(pow(x2, 2) + pow(y2, 2));

		if (d1 == 0 || d2 == 0)
		{
			continue;
		}

		double cosa = (x1*x2 + y1*y2) / (d1*d2);//这里应该是此向量和竖直向上的向量夹角的余弦值
		if (d >= 0)
		{
			NiCosS.push_back(cosa);
			NiIndexS.push_back(i);
			NiDis.push_back(d2);
		}
		else
		{
			CosS.push_back(cosa);
			IndexS.push_back(i);
			Dis.push_back(d2);	//存储距离
		}
	}
	if ((int)NiCosS.size() >= 1)
	{
		int Good;
		GetMinD(NiCosS, NiDis, Good);
		Good = NiIndexS[Good];
		Temppoint = m_point[Good];
		m_point[Good].mark = 1;
	}
	else
	{
		int Good;
		GetMinD(CosS, Dis, Good);
		Good = IndexS[Good];
		Temppoint = m_point[Good];
		m_point[Good].mark = 1;
	}
	Point.push_back(Temppoint);
	m_point[StarIndex].mark = 0;	//删除首点被使用痕迹
	for (int j = 2; j < Num + 1; j++)
	{
		NiIndexS.clear(); IndexS.clear();
		NiDis.clear();	  Dis.clear();
		NiCosS.clear();	  CosS.clear();
		for (int i = 0; i < Num; i++)
		{
			if (m_point[i].mark == 1)
			{
				continue;
			}
			double x2 = m_point[i].X - Point[j - 1].X;
			double y2 = m_point[i].Y - Point[j - 1].Y;
			double x1 = Point[j - 1].X - Point[j - 2].X;
			double y1 = Point[j - 1].Y - Point[j - 2].Y;
			double d = x1 * y2 - y1 * x2;		//叉乘结果
			double d1 = sqrt(pow(x1, 2) + pow(y1, 2));
			double d2 = sqrt(pow(x2, 2) + pow(y2, 2));
			if (d1 == 0 || d2 == 0)
			{
				continue;
			}
			double cosa = (-x1*x2 - y1*y2) / (d1*d2);
			if (d > 0)
			{
				NiCosS.push_back(cosa);
				NiIndexS.push_back(i);
				NiDis.push_back(d2);
			}
		}
		if ((int)NiCosS.size() >= 1)
		{
			int Good;
			GetMinD(NiCosS, NiDis, Good);
			Good = NiIndexS[Good];
			Temppoint = m_point[Good];
			m_point[Good].mark = 1;
		}
		else
		{
			return Point;
		}
		Point.push_back(Temppoint);
		if (Point[j].X == Point[0].X && Point[j].Y == Point[0].Y)
		{
			return Point;
		}
	}
	return Point;
}
vector<Data> Jie_Ds::LoadModel(const char* sFileName)
{
	vector<Data> Convexhull;
	std::ifstream in(sFileName);
	std::ifstream in2(sFileName);
	if (!in)
	{
		std::cout << "error of load mesh" << std::endl;
	}
	string str;
	Vector3D point2;

	in >> str;
	if (str == "v")
	{
		Data Points;
		in >> point2.x >> point2.y>>point2.z;
		//cout << point2.x << " " << point2.y << " "<<point2.z<<endl;
		Points.X = point2.x;
		Points.Y = point2.y;
		Points.Z = point2.z;
		Convexhull.push_back(Points);
	}
	while (getline(in, str))
		//while (!in.eof())
	{
		in >> str;
		if (str == "v")
		{
			Data Points;
			in >> point2.x >> point2.y>>point2.z;
			//cout << point2.x << " " << point2.y << " "<<point2.z<<endl;
			Points.X = point2.x;
			Points.Y = point2.y;
			Points.Z = point2.z;
			Convexhull.push_back(Points);
		}

	}
	return Convexhull;
}
void Jie_Ds::Charact()	
{
	double XMin = m_oriData[0].X;//m_oriData中存储原始数据，井位信息
	double XMax = m_oriData[0].X;
	double YMin = m_oriData[0].Y;
	double YMax = m_oriData[0].Y;
	double ZMin = m_oriData[0].Z;
	double ZMax = m_oriData[0].Z;
	for (size_t i = 1; i<m_oriData.size(); i++)
	{
		if (XMin > m_oriData[i].X)
			XMin = m_oriData[i].X;

		if (XMax < m_oriData[i].X)
			XMax = m_oriData[i].X;

		if (YMin > m_oriData[i].Y)
			YMin = m_oriData[i].Y;

		if (YMax < m_oriData[i].Y)
			YMax = m_oriData[i].Y;

		if (ZMin > m_oriData[i].Z)
			ZMin = m_oriData[i].Z;

		if (ZMax < m_oriData[i].Z)
			ZMax = m_oriData[i].Z;
	}

	/*************************************************/
	//是否要判断m_boder的数据？20131112  之前已经有这么一步的赋值操作：m_OriBoder = m_Border;
	for (size_t i = 0; i<m_Border.size(); i++)
	{
		if (XMin > m_Border[i].X)
			XMin = m_Border[i].X;

		if (XMax < m_Border[i].X)
			XMax = m_Border[i].X;

		if (YMin > m_Border[i].Y)
			YMin = m_Border[i].Y;

		if (YMax < m_Border[i].Y)
			YMax = m_Border[i].Y;
	}
	/*************************************************/

	m_XMin = XMin;
	m_XMax = XMax;
	m_YMin = YMin;
	m_YMax = YMax;
	m_ZMin = ZMin;
	m_ZMax = ZMax;
}
void Jie_Ds::DataOpt()
{
	vector<Data> datatemp;
	for (int i = 0; i < m_oriData.size(); i++)
	{
		Data p;
		p.X = m_oriData[i].X-m_XMin;
		p.Y = m_oriData[i].Y-m_YMin;
		p.Z = m_oriData[i].Z;
		datatemp.push_back(p);
	}
	m_oriData.clear();
	m_oriData = datatemp;
}
void Jie_Ds::DataRec()
{
	vector<Data> datatemp;
	for (int i = 0; i < m_oriData.size(); i++)
	{
		Data p;
		p.X = m_oriData[i].X + m_XMin;
		p.Y = m_oriData[i].Y + m_YMin;
		p.Z = m_oriData[i].Z;
		datatemp.push_back(p);
	}
	m_oriData.clear();
	m_oriData = datatemp;
}

void Jie_Ds::CalcBorder()
{
	if (m_Border.size() < 3) // 没有设置边界，计算一个边界
	{
		Data dt;
		dt.X = m_XMin; // 左上
		dt.Y = m_YMax;
		m_Border.push_back(dt);

		dt.X = m_XMax; // 右上
		dt.Y = m_YMax;
		m_Border.push_back(dt);

		dt.X = m_XMax; // 右下
		dt.Y = m_YMin;
		m_Border.push_back(dt);

		dt.X = m_XMin; // 左下
		dt.Y = m_YMin;
		m_Border.push_back(dt);
	}


	//m_Border.clear();//20131111
	//m_Border.resize(0);

	//for (int i = 0; i< m_Border.size(); i++)
	//{
	//	m_Border.push_back(m_Border[i]); // 导入边界
	//}

}
double Jie_Ds::GetDis(double x1, double y1, double x2, double y2)
{
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*((y1 - y2)));
}
void Jie_Ds::AddPt(vector<Data> &convexBag)
{
	int Num = (int)convexBag.size();
	double MaxD = 0;
	int Index = 0;
	for (int i = 0; i < Num - 1; i++)
	{
		Data Star = convexBag[i];
		Data End = convexBag[i + 1];
		double d = pow(Star.X - End.X, 2) + pow(Star.Y - End.Y, 2);
		if (d>MaxD)
		{
			MaxD = d;
			Index = i;
		}
	}
	Data Middle;
	Middle.X = 0.5*(convexBag[Index].X + convexBag[Index + 1].X);
	Middle.Y = 0.5*(convexBag[Index].Y + convexBag[Index + 1].Y);
	convexBag.insert(convexBag.begin() + Index + 1, Middle);
}
void Jie_Ds::deCasteljau()
{
	vector<point> bezierpoint;
	int n = control_point.size();

	for (double t = 0.0; t <= 1.0; t += 0.001 / n)//保留输入点
	{
		for (int i = 1; i < n; i++)
		{
			for (int j = 0; j < n - i; j++)
			{
				if (i == 1)//i=1由已知控制顶点计算
				{
					control_point[j].x = (1 - t) * input_vertice[j].x + t * input_vertice[j + 1].x;
					control_point[j].y = (1 - t) * input_vertice[j].y + t * input_vertice[j + 1].y;
					continue;
				}
				else//i!=1由上一次迭代结果计算
				{
					control_point[j].x = (1 - t) * control_point[j].x + t * control_point[j + 1].x;
					control_point[j].y = (1 - t) * control_point[j].y + t * control_point[j + 1].y;
				}
			}
		}
		b_spline.push_back(control_point[0]);
	}
}
void Jie_Ds::deBoor()
{
	float t[maxn];
	int k = 3;
	int n = control_point.size() - 1;
	//准均匀B样条
	for (int i = 0; i < n + k; i++)
	{
		if (i <= k - 1) t[i] = 0;
		if (i >= k && i < n + 1) t[i] = t[i - 1] + 1.0 / (n - k + 2);
		if (i >= n + 1) t[i] = 1;
	}
	for (int j = k - 1; j <= n; j++)
	{
		for (double u = t[j]; u <= t[j + 1]; u += 0.02 / n)
		{
			for (int r = 1; r <= k - 1; r++)
			{
				for (int i = j; i >= j - k + r + 1; i--)
				{
					float x1 = u - t[i];
					float x2 = t[i + k - r] - t[i];
					float y1 = t[i + k - r] - u;
					float k1, k2;
					if (x1 == 0.0 && x2 == 0.0) k1 = 0;
					else k1 = x1 / x2;
					if (y1 == 0.0 && x2 == 0.0) k2 = 0;
					else k2 = y1 / x2;

					if (r == 1)//控制点
					{
						control_point[i].x = input_vertice[i].x * k1 + input_vertice[i - 1].x * k2;
						control_point[i].y = input_vertice[i].y * k1 + input_vertice[i - 1].y * k2;
						continue;
					}
					else
					{
						control_point[i].x = control_point[i].x * k1 + control_point[i - 1].x * k2;
						control_point[i].y = control_point[i].y * k1 + control_point[i - 1].y * k2;
					}
				}
			}
			b_spline.push_back(control_point[j]);//递推的最后一层的点，即为求得的点
		}
	}


}
double Jie_Ds::bezier3funcX(double uu, Vector2D *controlP) {
	double part0 = controlP[0].x * uu * uu * uu;
	double part1 = 3 * controlP[1].x * uu * uu * (1 - uu);
	double part2 = 3 * controlP[2].x * uu * (1 - uu) * (1 - uu);
	double part3 = controlP[3].x * (1 - uu) * (1 - uu) * (1 - uu);
	return part0 + part1 + part2 + part3;
}
double Jie_Ds::bezier3funcY(double uu, Vector2D *controlP) {
	double part0 = controlP[0].y * uu * uu * uu;
	double part1 = 3 * controlP[1].y * uu * uu * (1 - uu);
	double part2 = 3 * controlP[2].y * uu * (1 - uu) * (1 - uu);
	double part3 = controlP[3].y * (1 - uu) * (1 - uu) * (1 - uu);
	return part0 + part1 + part2 + part3;
}
void Jie_Ds::createCurve()
{
	double scale = 0.6;
	int count = originPoint.size();
	//CvPoint midpoints[count];
	vector<Vector2D> midpoints;
	//生成中点       
	for (int i = 0; i < count; i++) {
		int nexti = (i + 1) % count;
		Vector2D p;
		p.x = (originPoint[i].x + originPoint[nexti].x) / 2.0;
		p.y = (originPoint[i].y + originPoint[nexti].y) / 2.0;
		midpoints.push_back(p);
	}

	//平移中点  
	Vector2D extrapoints[maxn];
	for (int i = 0; i < count; i++) {
		int nexti = (i + 1) % count;
		int backi = (i + count - 1) % count;
		Vector2D midinmid;
		midinmid.x = (midpoints[i].x + midpoints[backi].x) / 2.0;
		midinmid.y = (midpoints[i].y + midpoints[backi].y) / 2.0;
		double offsetx = originPoint[i].x - midinmid.x;
		double offsety = originPoint[i].y - midinmid.y;
		int extraindex = 2 * i;
		extrapoints[extraindex].x = midpoints[backi].x + offsetx;
		extrapoints[extraindex].y = midpoints[backi].y + offsety;
		//朝 originPoint[i]方向收缩   
		double addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;
		double addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;
		extrapoints[extraindex].x = originPoint[i].x + addx;
		extrapoints[extraindex].y = originPoint[i].y + addy;

		int extranexti = (extraindex + 1) % (2 * count);
		extrapoints[extranexti].x = midpoints[i].x + offsetx;
		extrapoints[extranexti].y = midpoints[i].y + offsety;
		//朝 originPoint[i]方向收缩   
		addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;
		addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;
		extrapoints[extranexti].x = originPoint[i].x + addx;
		extrapoints[extranexti].y = originPoint[i].y + addy;

	}

	Vector2D controlPoint[4];
	//生成4控制点，产生贝塞尔曲线  
	for (int i = 0; i < count; i++) {
		controlPoint[0] = originPoint[i];
		int extraindex = 2 * i;
		controlPoint[1] = extrapoints[extraindex + 1];
		int extranexti = (extraindex + 2) % (2 * count);
		controlPoint[2] = extrapoints[extranexti];
		int nexti = (i + 1) % count;
		controlPoint[3] = originPoint[nexti];
		double u = 1;
		while (u >= 0) {
			double px = bezier3funcX(u, controlPoint);
			double py = bezier3funcY(u, controlPoint);
			//u的步长决定曲线的疏密  
			u -= 0.01;
			Vector2D tempP = Vector2D(px, py);
			//存入曲线点   
			curvePoint.push_back(tempP);
		}
	}
}
void Jie_Ds::OptimizeBorderBezier(vector<Data> &convexBag,double e)
{
	std::ofstream outC("center_Bezier.obj");
	vector<Data> ConvexhullExp;
	double Sum_x = 0, Sum_y = 0, Aver_x = 0, Aver_y = 0;
	for (int j = 0; j < convexBag.size()-1; j++)
	{
		Sum_x += convexBag[j].X;
		Sum_y += convexBag[j].Y;
	}
	Aver_x = Sum_x / (convexBag.size() - 1);
	Aver_y = Sum_y / (convexBag.size() - 1);
	outC << "v " << fixed << setprecision(5) << Aver_x << " " << fixed << setprecision(5) << Aver_y << " " << "0" << endl;
	outC.close();
	for (int i = 0; i < convexBag.size(); i++)
	{
		Data cp;
		cp.X = convexBag[i].X;
		cp.Y = convexBag[i].Y;
		double d_v_x = convexBag[i].X - Aver_x;
		double d_v_y = convexBag[i].Y - Aver_y;
		double d_n_v_x = d_v_x / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		double d_n_v_y = d_v_y / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		cp.X += e* d_n_v_x;
		cp.Y += e* d_n_v_y;
		ConvexhullExp.push_back(cp);
	}

	for (int i = 0; i < ConvexhullExp.size()-1; i++)
	{
		Vector2D p;
		p.x = ConvexhullExp[i].X;
		p.y = ConvexhullExp[i].Y;
		originPoint.push_back(p);
	}
	createCurve();
	convexBag.clear();
	for (int i = 0; i < curvePoint.size(); i++)
	{
		Data d;
		d.X = curvePoint[i].x;
		d.Y = curvePoint[i].y;
		convexBag.push_back(d);
	}
	std::ofstream out101("convex_current_exp_Bezier.obj");
	for (int j = 0; j < ConvexhullExp.size(); j++)
	{
		Data qq;
		qq = ConvexhullExp[j];
		out101 << "v " << fixed << setprecision(5) << qq.X << " " << fixed << setprecision(5) << qq.Y << " " << "0" << endl;
	}
	int count_exp_ = 1;
	for (int j = 0; j < ConvexhullExp.size() - 1; j++)
	{
		out101 << "l " << count_exp_ << " " << count_exp_ + 1 << endl;
		count_exp_++;
	}
	out101.close();
	std::ofstream out2000("convex_bezier_Bezier.obj");
	int count_exp = 1;
	for (int j = 0; j < convexBag.size(); j++)
	{
		Data qq;
		qq = convexBag[j];
		out2000 << "v " << fixed << setprecision(5)<< qq.X << " " << fixed << setprecision(5)<<qq.Y << " " << "0" << endl;
	}
	count_exp_ = 1;
	for (int j = 0; j < convexBag.size() - 1; j++)
	{
		out2000 << "l " << count_exp << " " << count_exp + 1 << endl;
		count_exp++;
	}
	out2000 << "l " << convexBag.size() << " " << "1" << endl;
	out2000.close();


}
void Jie_Ds::OptimizeBoder(vector<Data> &convexBag,double e)
{
	vector<Data> ConvexhullExp;
	double Sum_x = 0, Sum_y = 0, Aver_x = 0, Aver_y = 0;
	for (int j = 0; j < convexBag.size()-1; j++)
	{
		Sum_x += convexBag[j].X;
		Sum_y += convexBag[j].Y;
	}
	Aver_x = Sum_x / (convexBag.size()-1);
	Aver_y = Sum_y / (convexBag.size()-1);
	for (int i = 0; i < convexBag.size(); i++)
	{
		Data cp;
		cp.X = convexBag[i].X;
		cp.Y = convexBag[i].Y;
		double d_v_x = convexBag[i].X - Aver_x;
		double d_v_y = convexBag[i].Y - Aver_y;
		double d_n_v_x = d_v_x / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		double d_n_v_y = d_v_y / sqrt((d_v_x*d_v_x) + (d_v_y*d_v_y));
		cp.X += e* d_n_v_x;
		cp.Y += e* d_n_v_y;
		ConvexhullExp.push_back(cp);
	}



	for (int i = 0; i < ConvexhullExp.size(); i++)
	{
		point p;
		p.x = ConvexhullExp[i].X;
		p.y = ConvexhullExp[i].Y;
		input_vertice.push_back(p);
	}
	control_point = input_vertice;
	deBoor();
	//deCasteljau();
	convexBag.clear();
	for (int i = 0; i < b_spline.size(); i++)
	{
		Data d;
		d.X = b_spline[i].x;
		d.Y = b_spline[i].y;
		convexBag.push_back(d);
	}
	std::ofstream out100("convex_current_exp.obj");
	for (int j = 0; j < ConvexhullExp.size(); j++)
	{
		Data qq;
		qq = ConvexhullExp[j];
		out100 << "v " << fixed << setprecision(5) << qq.X << " " << fixed << setprecision(5) << qq.Y << " " << "0" << endl;
	}
	int count_exp = 1;
	for (int j = 0; j < ConvexhullExp.size() - 1; j++)
	{
		out100 << "l " << count_exp << " " << count_exp + 1 << endl;
		count_exp++;
	}
	out100.close();



	std::ofstream out200("convex_bezier.obj");
	count_exp = 1;
	for (int j = 0; j < convexBag.size(); j++)
	{
		Data qq;
		qq = convexBag[j];
		out200 << "v " << qq.X << " " << qq.Y << " " << "0" << endl;
	}
	count_exp = 1;
	for (int j = 0; j < convexBag.size() - 1; j++)
	{
		out200 << "l " << count_exp << " " << count_exp + 1 << endl;
		count_exp++;
	}
	out200 << "l " << convexBag.size() << " " << "1" << endl;
	out200.close();

	//if ((int)convexBag.size() == 1)
	//{
	//	//制造四边形
	//	m_XMin = m_XMin - 200;
	//	m_XMax = m_XMax + 200;
	//	m_YMin = m_YMin - 200;
	//	m_YMax = m_YMax + 200;

	//	convexBag.clear();

	//	Data p1, p2, p3, p4;
	//	p1.X = m_XMin;
	//	p1.Y = m_YMin;

	//	p2.X = m_XMin;
	//	p2.Y = m_YMax;

	//	p3.X = m_XMax;
	//	p3.Y = m_YMax;

	//	p4.X = m_XMax;
	//	p4.Y = m_YMin;

	//	convexBag.push_back(p1);
	//	convexBag.push_back(p2);
	//	convexBag.push_back(p3);
	//	convexBag.push_back(p4);//左上角为（1,1）点，实际上的（0,0）点
	//	convexBag.push_back(p1);
	//}
	//else if ((int)convexBag.size() == 2)
	//{
	//	//制造对应的矩形
	//	convexBag.clear();

	//	Data p1, p2, p3, p4;
	//	if (m_XMin == m_XMax)
	//	{
	//		m_XMin = m_XMin - 200;
	//	}
	//	if (m_YMin == m_YMax)
	//	{
	//		m_YMin = m_YMin - 200;
	//	}
	//	p1.X = m_XMin;
	//	p1.Y = m_YMin;

	//	p2.X = m_XMin;
	//	p2.Y = m_YMax;

	//	p3.X = m_XMax;
	//	p3.Y = m_YMax;

	//	p4.X = m_XMax;
	//	p4.Y = m_YMin;

	//	convexBag.push_back(p1);
	//	convexBag.push_back(p2);
	//	convexBag.push_back(p3);
	//	convexBag.push_back(p4);		 //左上角为（1,1）点，实际上的（0,0）点
	//	convexBag.push_back(p1);        //这里为啥左上角有区别

	//}
	//else if ((int)convexBag.size() == 4)
	//{
	//	//添加一个点
	//	AddPt(convexBag);
	//}

	//double dx = m_XMax - m_XMin;
	//double dy = m_YMax - m_YMin;
	//if (dy - 8 * dx > 0)
	//{
	//	//修改图形为矩形？
	//}
	//else if (dx - 8 * dy > 0)
	//{
	//	//修改图形为矩形？
	//}
}
double Jie_Ds::Angle(Data &p0, const Data &p1, const Data &p2)
{
	//将第二个点沿着x轴方向向右平移一个单位赋值给第一个点
	p0.X = p1.X + 1;
	p0.Y = p1.Y;


	/*std::ofstream out("current_triangle_point.obj");
	out << "v " << p0.X << " " << p0.Y << " " << "0" << endl;
	out << "v " << p1.X << " " << p1.Y << " " << "0" << endl;
	out << "v " << p2.X << " " << p2.Y << " " << "0" << endl;

	out.close();*/

	double A, B, X, cross, angle;
	//计算 距离
	A = sqrt((p1.X - p0.X) * (p1.X - p0.X) + (p1.Y - p0.Y) *(p1.Y - p0.Y));
	B = sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) *(p1.Y - p2.Y));

	if (A == 0 || B == 0)
		return 0;

	//余弦值
	X = ((p0.X - p1.X) * (p2.X - p1.X) + (p0.Y - p1.Y) * (p2.Y - p1.Y)) / (A*B);

	//叉乘判断逆时针是否大于180度
	//a = (x1, y1) b = (x2, y2)
	//a×b = x1y2 - x2y1,若结果为正，则向量b在a的逆时针方向,否则，b在a的顺时针方向,若结果为0，则a与b共线
	cross = (p0.X - p1.X) * (p2.Y - p1.Y) - (p2.X - p1.X) * (p0.Y - p1.Y);

	if (X == 1 || X >1 || X < -1)
		return 0;

	double temp;
	double p = 3.14159265358979323846;
	if (X == -1)
	{
		angle = p;
	}
	else
	{
		if (cross < 0)//顺时针
		{
			//cout << "顺时针" << endl;
			//angle = atan2((p2.Y - p1.Y), (p2.X - p1.X));
			temp = (atan(-X / sqrt(-X * X + 1)) + 2 * atan(1.0));
			angle = p*2  - temp;
		}
		else//逆时针
		{
			//cout << "逆时针" << endl;
			//angle = atan2((p1.Y - p2.Y), (p2.X - p1.X));
			angle = (atan(-X / sqrt(-X * X + 1)) + 2 * atan(1.0));
		}
	}
	//if (angle > p)
	//	angle = angle - p;
	return angle;
}
void Jie_Ds::DividedFourParts(int n, vector<int>& FourPort)
{
	if (n>3)
	{
		FourPort[0] = n;
		FourPort[1] = (int)(n / 4);
		FourPort[2] = 2 * (int)(n / 4);
		FourPort[3] = n - (int)(n / 4);
	}
	else
	{
		FourPort[0] = 3;
		FourPort[1] = 1;
		FourPort[2] = 2;
		FourPort[3] = 2;
	}

}
void Jie_Ds::BordersPointDis(vector<int>& FourPart, vector<Data>& convexBag, vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
	vector<Dis>& upDis, vector<Dis>&downDis, vector<Dis>&leftDis, vector<Dis>&rightDis)
{
	//这个容量分配有什么依据
	downline.resize(FourPart[0] - FourPart[3] + 2);//
	leftLine.resize(FourPart[1] + 2);//
	upLine.resize(FourPart[2] - FourPart[1] + 2);//
	rightLine.resize(FourPart[3] - FourPart[2] + 2);//

	downDis.resize(downline.size());
	leftDis.resize(leftLine.size());
	upDis.resize(upLine.size());
	rightDis.resize(rightLine.size());

	bool JudgeMax = true;
	vector<int> side(4);
	//赋初值
	side[0] = FourPart[3];
	side[1] = 0;
	side[2] = FourPart[1];
	side[3] = FourPart[2];

	int count_s0 = 0, count_s1 = 0, count_s2 = 0, count_s3 = 0;
	int j = 1, i;
	while (JudgeMax)
	{
		i = 0;

		//给下边计数器小于最大值
		if (side[0] <= FourPart[0])
		{
			count_s0++;
			//给点赋值
			downline[j].X = convexBag[side[0]].X;
			downline[j].Y = convexBag[side[0]].Y;


			//当为首点时计算的距离 默认为到(0,0)的距离 则
			//为首点时 给距离赋值 0
			if (j == 1)
				downDis[j].dis = 0;
			else
			//计算此点距此边首点的距离
			downDis[j].dis = sqrt((downline[j].X - downline[j - 1].X) * (downline[j].X - downline[j - 1].X)
				+ (downline[j].Y - downline[j - 1].Y) * (downline[j].Y - downline[j - 1].Y)) + downDis[j - 1].dis;

			

			//边序号自增
			side[0] = side[0] + 1;
			//记录判断的数自增
			i = i + 1;
			/*if (side[0] == FourPart[0])
			{
				int j_ = j + 1;
				downline[j_].X = convexBag[0].X;
				downline[j_].Y = convexBag[0].Y;
				downDis[j_].dis = sqrt((downline[j_].X - downline[j_ - 1].X) * (downline[j_].X - downline[j_ - 1].X)
					+ (downline[j_].Y - downline[j_ - 1].Y) * (downline[j_].Y - downline[j_ - 1].Y)) + downDis[j_ - 1].dis;
			}*/

		}

		//左边赋值  同上
		if (side[1] <= FourPart[1])
		{
			count_s1++;
			leftLine[j].X = convexBag[side[1]].X;
			leftLine[j].Y = convexBag[side[1]].Y;

			if (j == 1)
				leftDis[j].dis = 0;
			else
			leftDis[j].dis = sqrt((leftLine[j].X - leftLine[j - 1].X) * (leftLine[j].X - leftLine[j - 1].X)
				+ (leftLine[j].Y - leftLine[j - 1].Y) * (leftLine[j].Y - leftLine[j - 1].Y)) + leftDis[j - 1].dis;

			

			side[1] = side[1] + 1;
			i = i + 1;
		}
		//上边赋值 同上
		if (side[2] <= FourPart[2])
		{
			count_s2++;
			upLine[j].X = convexBag[side[2]].X;
			upLine[j].Y = convexBag[side[2]].Y;
			if (j == 1)
				upDis[j].dis = 0;
			else
			upDis[j].dis = sqrt((upLine[j].X - upLine[j - 1].X) * (upLine[j].X - upLine[j - 1].X)
				+ (upLine[j].Y - upLine[j - 1].Y) * (upLine[j].Y - upLine[j - 1].Y)) + upDis[j - 1].dis;

			

			side[2] = side[2] + 1;
			i = i + 1;
		}

		//右边赋值 同上
		if (side[3] <= FourPart[3])
		{
			count_s3++;
			rightLine[j].X = convexBag[side[3]].X;
			rightLine[j].Y = convexBag[side[3]].Y;

			if (j == 1)
				rightDis[j].dis = 0;
			else
			rightDis[j].dis = sqrt((rightLine[j].X - rightLine[j - 1].X) * (rightLine[j].X - rightLine[j - 1].X)
				+ (rightLine[j].Y - rightLine[j - 1].Y) * (rightLine[j].Y - rightLine[j - 1].Y)) + rightDis[j - 1].dis;

			

			side[3] = side[3] + 1;
			i = i + 1;
		}

		//边的记录器 自增
		j = j + 1;
		//如果没有运算 i=0 则说明已经遍历完了 所有点上的点 跳出循环
		if (i == 0)
			JudgeMax = false;
	}
	std::ofstream outd("DP.obj");
	for (int i = 1; i < downline.size(); i++)
		outd << "v " << downline[i].X << " " << downline[i].Y << " " << "0" << endl;
	outd.close();
	std::ofstream outu("UP.obj");
	for (int i = 1; i < upLine.size(); i++)
		outu << "v " << upLine[i].X << " " << upLine[i].Y << " " << "0" << endl;
	outu.close();
	std::ofstream outl("LP.obj");
	for (int i = 1; i < leftLine.size(); i++)
		outl << "v " << leftLine[i].X << " " << leftLine[i].Y << " " << "0" << endl;
	outl.close();
	std::ofstream outr("RP.obj");
	for (int i = 1; i < rightLine.size(); i++)
		outr << "v " << rightLine[i].X << " " << rightLine[i].Y << " " << "0" << endl;
	outr.close();
	//cout << "j = :" << j << endl;
	//cout << "count_s0 = :" << count_s0 << endl;
	//cout << "count_s1 = :" << count_s1 << endl;
	//cout << "count_s2 = :" << count_s2 << endl;
	//cout << "count_s3 = :" << count_s3 << endl;
}
//计算upLine、downline、leftLine、rightLine中Dis分量占总的长度的半分比
//计算upDis、downdis、leftDis、rightDis中angle分量，和前一个点构成的线段与x轴正方向的夹角
void Jie_Ds::BordersChar(vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
	vector<Dis>& upDis, vector<Dis>&downdis, vector<Dis>&leftDis, vector<Dis>&rightDis)
{
	bool JudgeMax = true;
	Data triPoint[3];
	int i, j = 2;
	std::ofstream out22("2.obj");
	while (JudgeMax)
	{
		//一个判断次循环是否有计算的 计数器 i>0说明次循环有运算，否则i=0 跳出循环
		//赋初值
		i = 0;
		//个数没有超出 下边 则进行计算
		if (j < (int)downline.size())
		{
			//计算此点占此边总长的百分数
			//这里downline.size()比downdis.size()大1，downsize下标从1开始，downsize[1]=0;
			downdis[j].per = downdis[j].dis / downdis[downline.size() - 1].dis;
			//赋值记录此点与前一点的线段
			triPoint[1] = downline[j - 1];
			triPoint[2] = downline[j];

			//计算此线段与x轴正方向的角度
			downdis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//cout << downdis[j].angle << endl;
			//自增 >0说明有运算
			i = i + 1;
			
			out22 << "v " << downline[downline.size() - 1].X << " " << downline[downline.size() - 1].Y << " " << "0" << endl;
		}
		
		if (j < (int)leftLine.size())
		{
			//计算此点占此边总长的百分数
			leftDis[j].per = leftDis[j].dis / leftDis[leftLine.size() - 1].dis;
			//赋值记录此点与前一点的线段
			triPoint[1] = leftLine[j - 1];
			triPoint[2] = leftLine[j];

			//计算此线段与x轴正方向的角度
			leftDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//自增 >0说明有运算
			i = i + 1;
			
		}

		if (j < (int)upLine.size())
		{
			//计算此点占此边总长的百分数
			upDis[j].per = upDis[j].dis / upDis[upLine.size() - 1].dis;
			//赋值记录此点与前一点的线段
			triPoint[1] = upLine[j - 1];
			triPoint[2] = upLine[j];

			//计算此线段与x轴正方向的角度
			upDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//自增 >0说明有运算
			i = i + 1;
			
		}

		if (j < (int)rightLine.size())
		{
			//计算此点占此边总长的百分数
			rightDis[j].per = rightDis[j].dis / rightDis[rightLine.size() - 1].dis;
			//赋值记录此点与前一点的线段
			triPoint[1] = rightLine[j - 1];
			triPoint[2] = rightLine[j];

			//计算此线段与x轴正方向的角度
			rightDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//自增 >0说明有运算
			i = i + 1;
			
		}

		//边界计数器自增
		j = j + 1;
		//i=0则说明此次循环没有运算 则跳出循环
		if (i == 0)
			JudgeMax = false;
		out22.close();

	}
}
void Jie_Ds::BordersPoints(vector<int>& FourPart, vector<Data>& convexBag, vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
	vector<Dis>& upDis, vector<Dis>&downdis, vector<Dis>&leftDis, vector<Dis>&rightDis,
	vector<Data>& DL, vector<Data>& DR, vector<Data>&DD, vector<Data>& DU,
	vector<TWOVALUE>& sim1, vector<TWOVALUE>& sim2, vector<TWOVALUE>& sim3, vector<TWOVALUE>& sim4, int XM, int YM)
{
	int i, j, m;
	double s, t, Tran;

	//*****************分割算法
	//确定四点坐标
	//存储映射的原点坐标
	sim1[1].X = convexBag[FourPart[0]].X;//这里应该存储的是凸包边界上最后一个点的位置

	//对应第一条边的x坐标差值
	sim1[2].X = convexBag[FourPart[3]].X - convexBag[FourPart[0]].X;

	//对应第二条边的x坐标差值
	sim1[3].X = convexBag[FourPart[1]].X - convexBag[FourPart[0]].X;

	sim1[4].X = convexBag[FourPart[0]].X - convexBag[FourPart[1]].X + convexBag[FourPart[2]].X - convexBag[FourPart[3]].X;

	//解释同上 对应y坐标值
	sim1[1].Y = convexBag[FourPart[0]].Y;
	sim1[2].Y = convexBag[FourPart[3]].Y - convexBag[FourPart[0]].Y;
	sim1[3].Y = convexBag[FourPart[1]].Y - convexBag[FourPart[0]].Y;
	sim1[4].Y = convexBag[FourPart[0]].Y - convexBag[FourPart[1]].Y + convexBag[FourPart[2]].Y - convexBag[FourPart[3]].Y;

	/*cout << "=======================" << endl;
	cout << "=======================" << endl;
	for(int p=1;p<=4;p++)
	cout << "v " << sim1[p].X << " " << sim1[p].Y << " " << "0" << endl;*/



	sim3.resize(YM + 1);
	sim2.resize(XM + 1);
	sim4.resize(XM + 1);
	DL.resize(YM + 1);
	DR.resize(YM + 1);
	DD.resize(XM + 1);
	DU.resize(XM + 1);
	for (j = 1; j <= YM; j++)
	{
		//此分割点距边起点的距离占此边长度的百分数
		t = double(j - 1) / (YM - 1);//写的应该有问题吧？？？？？

									 //计算
									 //则此点到此边起点的距离
		sim3[j].X = sim1[3].X * t;
		sim3[j].Y = sim1[3].Y * t;

		//遍历此点的特征数组(存储 点到边的距离,及此距离为此边的百分数, _和此点与前一点的线段和x轴正方向的角度)
		for (m = 2; m<(int)leftDis.size(); m++)
		{
			//当分割点 百分数与 凸包数组的百分数相同时
			if (t == leftDis[m - 1].per) //当凸包点为分割点时
			{
				//把此凸包点赋值给左边分割点数组
				DL[j].X = leftLine[m - 1].X;
				DL[j].Y = leftLine[m - 1].Y;
				//查找到此分割点遍历结束
				break;
			}

			//当分割点遍历到此边的尾点时
			else if (t == 1)
			{
				//把此边的最后一个点赋值给 分割点数组
				DL[j].X = convexBag[FourPart[1]].X;
				DL[j].Y = convexBag[FourPart[1]].Y;
				//查找到了跳出循环
				break;
			}

			//当此分割点距首点的距离占的百分数 为在条凸包边上时 则  t在某个线段两端点的百分数之间
			else if (leftDis[m - 1].per < t && t < leftDis[m].per)
			{
				//此分割点所在凸包边上的位置到此凸包边一端点距离占此凸包边长度的百分数
				Tran = leftDis[leftDis.size() - 1].dis* t - leftDis[m - 1].dis;
				//根据此边的角度与此点到一段点的距离 由直线参数方程 可计算此点的(x,y)坐标
				DL[j].X = leftLine[m - 1].X + Tran * cos(leftDis[m].angle);
				DL[j].Y = leftLine[m - 1].Y + Tran * sin(leftDis[m].angle);
				break;
			}
		}

		if ((int)rightDis.size()>1)
		{
			//同上 遍历右边界的分割点 并给它赋值
			for (m = 2; m< (int)rightDis.size(); m++)
			{
				if ((1 - t) == rightDis[m - 1].per)  //当为某个节点时
				{
					DR[j].X = rightLine[m - 1].X;
					DR[j].Y = rightLine[m - 1].Y;
					break;
				}
				else if ((1 - t) == 1)
				{
					DR[j].X = convexBag[FourPart[3]].X;
					DR[j].Y = convexBag[FourPart[3]].Y;
					break;
				}
				else if (rightDis[m - 1].per < (1 - t) && (1 - t) < rightDis[m].per)
				{
					Tran = rightDis[rightDis.size() - 1].dis * (1 - t) - rightDis[m - 1].dis;
					DR[j].X = rightLine[m - 1].X + Tran * cos(rightDis[m].angle);
					DR[j].Y = rightLine[m - 1].Y + Tran * sin(rightDis[m].angle);
					break;
				}
			}
		}
		//当特征数组=3时, 只有一个点
		else
		{
			//此点为第4个点
			DR[j].X = convexBag[FourPart[3]].X;
			DR[j].Y = convexBag[FourPart[3]].Y;
		}
	}

	// 左右边界上分割点遍历完毕

	//遍历上下边界的分割点
	for (i = 1; i <= XM; i++)
	{
		//计算此点占总边长的百分数
		s = double(i - 1) / (XM - 1);
		sim2[i].X = sim1[2].X * s;
		sim2[i].Y = sim1[2].Y * s;
		sim4[i].X = sim1[4].X * s;
		sim4[i].Y = sim1[4].Y * s;

		//给上边界分割点赋值
		for (m = 2; m< (int)upDis.size(); m++)
		{
			//当此分割点的百分数与凸包点的百分数相同时 凸包点为分割点
			if (s == upDis[m - 1].per) //当为某个节点时
			{
				//给分割点赋值
				DU[i].X = upLine[m - 1].X;
				DU[i].Y = upLine[m - 1].Y;
				//成立时跳出循环
				break;
			}

			//为尾点时
			else if (s == 1)
			{
				DU[i].X = convexBag[FourPart[2]].X;
				DU[i].Y = convexBag[FourPart[2]].Y;
				//跳出循环
				break;

			}

			//当为某个凸包边时
			else if (upDis[m - 1].per < s && s < upDis[m].per)
			{
				Tran = upDis[upDis.size() - 1].dis * s - upDis[m - 1].dis;
				DU[i].X = upLine[m - 1].X + Tran * cos(upDis[m].angle);
				DU[i].Y = upLine[m - 1].Y + Tran * sin(upDis[m].angle);
				break;
			}
		}

		//与上同理 可计算下边界分割点的坐标值
		for (m = 2; m < (int)downdis.size(); m++)
		{
			if ((1 - s) == downdis[m - 1].per)  //当为某个节点时
			{
				DD[i].X = downline[m - 1].X;
				DD[i].Y = downline[m - 1].Y;
				break;
			}

			else if ((1 - s) == 1)
			{
				DD[i].X = convexBag[FourPart[0]].X;
				DD[i].Y = convexBag[FourPart[0]].Y;
				break;
			}
			else if (downdis[m - 1].per <= 1 - s && 1 - s <= downdis[m].per)
			{
				Tran = downdis[downdis.size() - 1].dis * (1 - s) - downdis[m - 1].dis;
				DD[i].X = downline[m - 1].X + Tran * cos(downdis[m].angle);
				DD[i].Y = downline[m - 1].Y + Tran * sin(downdis[m].angle);
				break;
			}

		}

	}
	std::ofstream outDL("FGL.obj");
	for (int i = 1; i < DL.size(); i++)
		outDL << "v " << DL[i].X << " " << DL[i].Y << " " << "0" << endl;
	outDL.close();
	cout << DR.size() << endl;
	std::ofstream outDR("FGR.obj");
	for (int i = 1; i < DR.size(); i++)
		outDR << "v " << DR[i].X << " " << DR[i].Y << " " << "0" << endl;
	outDR.close();
	std::ofstream outDU("FGU.obj");
	for (int i = 1; i < DU.size(); i++)
		outDU << "v " << DU[i].X << " " << DU[i].Y << " " << "0" << endl;
	outDU.close();
	std::ofstream outDD("FGD.obj");
	for (int i = 1; i < DD.size(); i++)
		outDD << "v " << DD[i].X << " " << DD[i].Y << " " << "0" << endl;
	outDD.close();


}
void Jie_Ds::AddData(Data &t)
{
	bool m_b = false;
	for (int i = 0; i<(int)m_oriData.size(); i++)
	{
		if (m_oriData[i].X == t.X && m_oriData[i].Y == t.Y)
			m_b = true;
	}
	if (!m_b)
		m_oriData.push_back(t);
}
bool Jie_Ds::Inv(vector<vector<double>>&M)
{
	int i, j, k, n;
	double Temp;
	n = (int)M.size();
	vector<int> iw(n), jw(n);
	//在m右边增加一个单位阵，构成一个m的增广矩阵mm
	//double **mm = new double * [n];
	//vector<vector<double>> mm(n);
	//for(i = 0 ;i< n;i++)
	//{
	//	mm[i].resize(2 * n);
	//	for(int j = 0;j < n;j ++) 
	//	{
	//		mm[i][j] = M[i][j];
	//	}
	//}
	double **mm = new double *[n];		//数据比较大或操作比较多，vector效率比指针低很多20150325
	for (i = 0; i< n; i++)
	{
		mm[i] = new double[2 * n];
		for (int j = 0; j < n; j++)
		{
			mm[i][j] = M[i][j];
		}
	}

	for (i = 0; i< n; i++)
	{
		for (j = n; j<2 * n; j++)
		{
			if (i == j - n)
				mm[i][j] = 1;
			else
				mm[i][j] = 0;
		}
	}
	//通过初等行变换(即高斯消去法)使原矩阵变为单位阵，则右边的单位阵即是原矩阵的逆阵
	for (k = 0; k < n - 1; k++)
	{
		/*----------------------------------------*/
		//从第 k 行、第 k 列开始的右下角子阵中选取绝对值最大的元素，并记住次元素在的行号和列号，
		//在通过行交换和列交换将它交换到主元素位置上.这一步称为全选主元  20140925
		iw[k] = k;
		jw[k] = k;
		if (abs(mm[k][k]) < 0.00000001)
		{
			for (int i = k + 1; i < n; i++)
			{
				if (abs(mm[i][k]) > 0.000001)
				{
					iw[k] = i;
					for (j = 0; j < n; j++)
					{
						Temp = mm[k][j];
						mm[k][j] = mm[iw[k]][j];
						mm[iw[k]][j] = Temp;
						//swap(mm[k][j], mm[iw[k]][j]);
					}
					break;
				}
			}
		}
		if (abs(mm[k][k]) < 0.00000001)
		{
			return false;
		}
		/*----------------------------------------*/
		for (i = k + 1; i < n; i++)
		{
			Temp = mm[i][k] / mm[k][k];
			for (j = 0; j < 2 * n; j++)
			{
				mm[i][j] = mm[i][j] - mm[k][j] * Temp;
			}
			mm[i][k] = 0.0;		//防止计算误差20140929
		}
	}

	if (abs(mm[n - 1][n - 1]) < 0.00000001)
	{
		return false;
	}
	for (k = n - 1; k > 0; k--)
	{
		for (i = k - 1; i >= 0; i--)
		{
			Temp = mm[i][k] / mm[k][k];
			for (j = 2 * n - 1; j >= 0; j--)
				mm[i][j] = mm[i][j] - mm[k][j] * Temp;
		}
	}
	double s;
	for (i = 0; i<n; i++)
	{
		s = mm[i][i];
		for (j = 0; j < 2 * n; j++)
			mm[i][j] = mm[i][j] / s;
	}

	//输出变换后的右边的矩阵
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			M[i][j] = mm[i][j + n];
	}

	//根据在全选主元过程中所记录的行、列交换的信息进行恢复，恢复的原则如下：在全选主元过程中，
	//先交换的行（列）后进行恢复；原来的行（列）交换用列（行）交换来恢复。
	for (k = n - 2; k >= 0; k--)
	{
		if (jw[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				Temp = M[k][i];
				M[k][i] = M[jw[k]][i];
				M[jw[k]][i] = Temp;
				//swap(MInv[k][i], MInv[jw[k]][i]);
			}
		}

		if (iw[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				Temp = M[i][k];
				M[i][k] = M[i][iw[k]];
				M[i][iw[k]] = Temp;
				//swap(MInv[i][k], MInv[i][iw[k]]);
			}
		}
	}
	return true;
}
void Jie_Ds::PreMatrix(vector<Data>& suV, vector<double>& ni)
{
	double b0 = 1;			//先假设斜率为一个值
	double b1 = 0;			//设截距为一个值
	double RR = m_B;		//设最大滞后距
	int n = (int)suV.size();
	vector<vector<double>> Arr(n + 1);

	//********给矩阵赋值
	//&&&&&&&&
	Arr[0].resize(n + 1);
	for (int i = 1; i <= n; i++)
	{
		Arr[i].resize(n + 1);
		Arr[i][i] = b1 + b0 * RR - b1;//0?
		Arr[i][n] = 1;
		Arr[0][i] = 1;
		for (int j = 0; j <= n - 1; j++)
		{
			double d = sqrt((suV[i - 1].X - suV[j].X) * (suV[i - 1].X - suV[j].X) + (suV[i - 1].Y - suV[j].Y) * (suV[i - 1].Y - suV[j].Y));
			if (d > RR)
			{
				d = RR;			//20140422最大滞后距离
			}
			Arr[i][j] = abs(b0 * (d)-RR);
		}
	}

	Arr[0][n] = 0;
	Arr[0][0] = 1;

	/************************************************/
	//CString Str = "z值 \n";
	//for (int i = 1;i <= n;i++)
	//{
	//	for(int j = 1;j <= n;j++)
	//	{
	//		CString s;
	//		s.Format("    %f",Arr[i][j]);
	//		Str += s;
	//	}
	//	Str += "\n";
	//}
	//for (int i = 0 ; i < suV.size() ; i ++)
	//{
	//	CString str;
	//	str.Format("%f,    %f,    %f，   %f     ，%f\n",suV[i].X,suV[i].Y,suV[i].Z,m_ZMax,m_ZMin);
	//	Str += str;
	//}
	//AfxMessageBox(Str);
	/************************************************/
	//求出矩阵的逆矩阵
	Inv(Arr);//求出矩阵的逆矩阵要是使用eigen库的话，只需要一行代码

			 //求解方程根
	ni.resize(n + 1);

	Data d;
	d.X = 0; d.Y = 0; d.Z = 0;
	suV.push_back(d);
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			ni[i] = ni[i] + suV[j].Z * Arr[j][i];
		}
	}
}
double Jie_Ds::InsertEst(vector<Data>& suV, TWOVALUE& D, vector<double>& ni)
{
	double dg = 0;
	double ZZ = m_B;//直线型变异函数B值
	double b0 = 1;
	double b1 = 0;
	vector<double> cc(suV.size());
	cc[0] = 1;
	dg = cc[0] * ni[0];
	double distanceVal = 0;
	for (int i = 1; i < (int)suV.size(); i++)
	{
		distanceVal = 0;
		distanceVal = sqrt((suV[i - 1].X - D.X) * (suV[i - 1].X - D.X) + (suV[i - 1].Y - D.Y) * (suV[i - 1].Y - D.Y));
		if (distanceVal > ZZ)
		{
			distanceVal = ZZ;	//20140422最大滞后距离
		}
		cc[i] = b0 * (ZZ - distanceVal);
		dg = cc[i] * ni[i] + dg;
	}
	//cout << dg << " " << endl;
	return dg;
}
double Jie_Ds::DisInv(TWOVALUE D)
{
	double Z = 0;
	double Sum = 0;
	double Di = 0;
	Data TempWell;
	for (int i = 0; i < (int)m_suV.size(); i++)
	{
		TempWell = m_suV[i];
		double d = (pow(D.X - TempWell.X, 2) + pow(D.Y - TempWell.Y, 2));
		//d = pow(d,1.5);		//加大距离反比法的次数，以便离得近的点影响力更大20131217
		if (d < 0.00001)
		{
			return TempWell.Z;
		}
		/************************/
		//TWOVALUE temp;
		//temp.X = TempWell.X;
		//temp.Y = TempWell.Y;
		//if (IsChangeW(D,temp,m_vecFaultLines))
		//{
		//	d = d * 1000;								//隔着断层则缩小权值20140904；
		//}
		/**************************/
		//d = pow(d,3) ;
		d = 1.0 / d;
		Sum += d;
		Di = Di + d * TempWell.Z;
	}
	Z = Di / Sum;
	return Z;
}
void Jie_Ds::EvaluateNoFault()
{
	if ((int)m_oriData.size() < 800)//20150323,800口井可以使用克里金方法，速度4s左右可以接受
	{
		m_IsK = true;				//点数较少使用克里金估值
	}
	else
	{
		m_IsK = false;				//点数太多使用幂距离反比估值
	}

	if ((int)m_oriData.size() <= 1)
	{
		return;		//处理只有一个井位点的情况20131107
	}
	if (m_IsK)			//是否使用克里金算法估值（准确但是速度慢）201308163
	{
		vector<Data> Suv = m_oriData;
		vector<double> ni;
		PreMatrix(Suv, ni);
		TWOVALUE D;
		for (int i = 1; i <= m_XNum; i++)
		{
			for (int j = 1; j <= m_YNum; j++)
			{
				D.X = m_GridPoint[i][j].X;
				D.Y = m_GridPoint[i][j].Y;
				//估计此点的属性值    没有断层的估计值 + 修正值
				m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);//全部井位点的克里金估值

															//m_GridPoint[i][j].Z = Well_Near_K(D);		//邻井范围内的克里金估值20131212
															//m_GridPoint[i][j].Z = DisInv(D);			//邻井范围内的距离反比估值20131212
				if (m_ZMin > m_GridPoint[i][j].Z)
					m_ZMin = m_GridPoint[i][j].Z;

				if (m_ZMax < m_GridPoint[i][j].Z)
					m_ZMax = m_GridPoint[i][j].Z;
			}
		}
		m_suV = Suv;		//20131212
		m_ni = ni;			//20130913 应用于外推法
	}
	else
	{
		m_suV = m_oriData;
		TWOVALUE D;
		for (int i = 1; i <= m_XNum; i++)
		{
			for (int j = 1; j <= m_YNum; j++)
			{
				D.X = m_GridPoint[i][j].X;
				D.Y = m_GridPoint[i][j].Y;
				//估计此点的属性值    没有断层的估计值 + 修正值
				//m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);//全部井位点的克里金估值
				//m_GridPoint[i][j].Z = Well_Near_K(D);		//邻井范围内的克里金估值20131212
				m_GridPoint[i][j].Z = DisInv(D);			//邻井范围内的距离反比估值20131212
				if (m_ZMin > m_GridPoint[i][j].Z)
					m_ZMin = m_GridPoint[i][j].Z;

				if (m_ZMax < m_GridPoint[i][j].Z)
					m_ZMax = m_GridPoint[i][j].Z;
			}

		}
	}

	//if ((int)m_oriData.size() <= 1)
	//{
	//	return ;		//处理只有一个井位点的情况20131107
	//}
	//vector<double> ni;
	//vector<Data> Suv = m_oriData; //20140421
	//PreMatrix(Suv,ni);
	//TWOVALUE D;
	//for (int i = 1;i <= m_XNum;i++)
	//{
	//	for(int j = 1;j <= m_YNum;j++)
	//	{
	//		D.X= m_GridPoint[i][j].X;
	//		D.Y = m_GridPoint[i][j].Y;
	//		//估计此点的属性值    没有断层的估计值 + 修正值
	//		m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);
	//		if(m_ZMin > m_GridPoint[i][j].Z)
	//			m_ZMin = m_GridPoint[i][j].Z;

	//		if(m_ZMax < m_GridPoint[i][j].Z)
	//			m_ZMax = m_GridPoint[i][j].Z;
	//	}

	//}

	/************************************************/
	//CString Str = "z值 \n";
	//for (int i = 1;i <= m_XNum;i++)
	//{
	//	for(int j = 1;j <= m_YNum;j++)
	//	{
	//		CString s;
	//		s.Format("    %f",m_GridPoint[i][j].Z);
	//		Str += s;
	//	}
	//	Str += "\n";
	//}
	//for (int i = 0 ; i < m_oriData.size() ; i ++)
	//{
	//	CString str;
	//	str.Format("%f,    %f,    %f，   %f     ，%f\n",m_oriData[i].X,m_oriData[i].Y,m_oriData[i].Z,m_ZMax,m_ZMin);
	//	Str += str;
	//}
	//AfxMessageBox(Str);
	/************************************************/
	//m_suV = Suv;
	//m_ni = ni;				//20130913 应用于外推法
}
void Jie_Ds::SetGridXY()
{
	//这里的m_Border中应该存储的是矩形边界
	m_OriBoder = m_Border;	//原始边界信息由m_OriBoder  20140806
	vector<Data> convexBag;
	//Charact();				//求取最大最小值

							/*-------------------------------------------*/
							//兼容二级边界不闭合的情况20170718
	int OriCount = (int)m_OriBoder.size();
	if (OriCount >= 3)
	{
		Data Pt0 = m_OriBoder[0];
		Data Pt1 = m_OriBoder[OriCount - 1];
		if (GetDis(Pt0.X, Pt0.Y, Pt1.X, Pt1.Y) >= 0.0001)
		{
			//让边界闭合
			m_OriBoder.push_back(Pt0);
		}
	}

	/*-------------------------------------------*/

	/*-------------------------------------------*/
	//根据井位距离计算变异距离m_B的值 20150205  变异距离是什么？？
	double k = sqrt(pow(m_XMax - m_XMin, 2) + pow(m_YMax - m_YMin, 2));
	if (m_B <= k)
	{
		m_B = k + 2;		//主要是防止西北那边特别大的范围
	}
	/*-------------------------------------------*/

	//GetRectBoder();		//使用矩形边界（自己构造边界）20131025

	//判断是否有边界；如果有边界，就从边界计算凸包，否则从数据点计算凸包
	//20131111m_Border中一直有数值，区别是单纯的井位数据还是人工数据
	//由数据点进行凸包的计算
	if ((int)m_Border.size() <= 2)
	{
		m_Border.clear();
		//convexBag = Convex(m_oriData);
		convexBag = Withershins(m_oriData);	//测试新的凸包函数20131108 构造凸包
		//使用B样条进行凸包边界的优化
		//OptimizeBoder(convexBag,200);			
		//使用贝塞尔曲线进行凸包边界的优化
		OptimizeBorderBezier(convexBag,1000);
		m_Border = convexBag;
		m_OriBoder = m_Border;	//原始边界信息由m_OriBoder  20140806
	}
	//由边界进行凸包的计算
	else
	{
		//cout << "hello2" << endl;
		convexBag = Withershins(m_Border);  //测试20131107  构造凸包
		//使用B样条进行凸包边界的优化
		//OptimizeBoder(convexBag);
		//使用贝塞尔曲线进行凸包边界的优化
		OptimizeBorderBezier(convexBag,3); 
		m_Border = convexBag;
	}

	std::ofstream out("convexhull_test_point.obj");
	//out << "v " << convexBag[379].X << " " << convexBag[376].Y << " " << "0" << endl;
	//out << "v " << convexBag[499].X << " " << convexBag[499].Y << " " << "0" << endl;
	/*out << "v " << convexBag[377].X << " " << convexBag[377].Y << " " << "0" << endl;
	out << "v " << convexBag[378].X << " " << convexBag[378].Y << " " << "0" << endl;
	out << "v " << convexBag[389].X << " " << convexBag[389].Y << " " << "0" << endl;
	out << "v " << convexBag[390].X << " " << convexBag[390].Y << " " << "0" << endl;
	out << "v " << convexBag[391].X << " " << convexBag[391].Y << " " << "0" << endl;*/
	//out << "v " << convexBag[0].X << " " << convexBag[0].Y << " " << "0" << endl;
	//out << "v " << convexBag[126].X << " " << convexBag[126].Y << " " << "0" << endl;
	//out << "v " << convexBag[252].X << " " << convexBag[252].Y << " " << "0" << endl;
	//out.close();
	vector<int> FourPort(4);
	

	DividedFourParts((int)(convexBag.size() - 1), FourPort);//优化边界后凸包点分为上下左右四组

	//cout << "=========================================" << endl;
	//cout << "输出凸包点分组的一些信息" << endl;
	//cout << "FourPort[0]= " << FourPort[0] << " " << "FourPort[1]= " << FourPort[1] << " " << "FourPort[2]= " << FourPort[2] << " " << "FourPort[3]= " << FourPort[3] << endl;

	vector<Data> upLine, downLine, leftLine, rightLine;

	vector<Dis> upDis, downdis, leftDis, rightDis;

	//upLine、downline、leftLine、rightLine存储了凸包四部分边界上的点的位置
	//upDis、downDis、leftDis、rightDis中的dis分量存储了对应点到第一个点的距离
	BordersPointDis(FourPort, convexBag, upLine, downLine, leftLine, rightLine,
		upDis, downdis, leftDis, rightDis);

	//计算upLine、downline、leftLine、rightLine中Dis分量占总的长度的半分比
	//计算upDis、downdis、leftDis、rightDis中angle分量，与x轴正方向的夹角
	BordersChar(upLine, downLine, leftLine, rightLine, upDis, downdis, leftDis, rightDis);

	vector<Data> DL, DR, DD, DU;
	vector<TWOVALUE> sim1(5), sim2, sim3, sim4;
	BordersPoints(FourPort, convexBag, upLine, downLine, leftLine, rightLine,
		upDis, downdis, leftDis, rightDis, DL, DR, DD, DU,
		sim1, sim2, sim3, sim4, m_XNum, m_YNum);

	double s, t;
	int XM = m_XNum, YM = m_YNum;

	int GridPoint_size = m_GridPoint.size();
	if (GridPoint_size >0)
	{
		for (int i = 0; i<GridPoint_size; i++)
		{
			vector <THRVALUE>().swap(m_GridPoint[i]);
			m_GridPoint[i].clear();
		}
	}
	vector <vector<THRVALUE>>().swap(m_GridPoint);
	m_GridPoint.clear();
	m_GridPoint.resize(m_XNum + 3);
	m_GridPoint[0].resize(m_XNum + 3);

	for (int i = 1; i <= m_XNum; i++)
	{
		m_GridPoint[i].resize(m_XNum + 3);
	}

	m_GridPoint[m_XNum + 1].resize(m_XNum + 3);
	m_GridPoint[m_XNum + 2].resize(m_XNum + 3);

	for (int i = 1; i <= XM; i++)
	{
		s = (double)(i - 1) / (XM - 1);
		for (int j = 1; j <= YM; j++)
		{
			t = (double)(j - 1) / (YM - 1);
			m_GridPoint[i][j].X = (DL[j].X * (1 - s) + DR[j].X * s) + (DU[i].X * t + DD[i].X * (1 - t)) - (sim1[1].X + sim2[i].X + sim3[j].X + t * sim4[i].X);
			m_GridPoint[i][j].Y = (DL[j].Y * (1 - s) + DR[j].Y * s) + (DU[i].Y * t + DD[i].Y * (1 - t)) - (sim1[1].Y + sim2[i].Y + sim3[j].Y + t * sim4[i].Y);
			m_GridPoint[i][j].Z = 0.0;//20130814
		}
	}

	DL.clear();//清除后再交换？
	vector <Data>().swap(DL);

	DR.clear();
	vector <Data>().swap(DR);

	DD.clear();
	vector <Data>().swap(DD);

	DU.clear();
	vector <Data>().swap(DU);

	sim2.clear();
	vector <TWOVALUE>().swap(sim2);

	sim3.clear();
	vector <TWOVALUE>().swap(sim3);

	sim4.clear();
	vector <TWOVALUE>().swap(sim4);

	sim1.clear();
	vector <TWOVALUE>().swap(sim1);

	upLine.clear();
	vector <Data>().swap(upLine);

	downLine.clear();
	vector <Data>().swap(downLine);

	leftLine.clear();
	vector <Data>().swap(leftLine);

	rightLine.clear();
	vector <Data>().swap(rightLine);

	upDis.clear();
	vector <Dis>().swap(upDis);

	downdis.clear();
	vector <Dis>().swap(downdis);

	leftDis.clear();
	vector <Dis>().swap(leftDis);

	rightDis.clear();
	vector <Dis>().swap(rightDis);
}
//得到一个数的数量级
double Jie_Ds::GetMagnitude(double fNumber)
{
	//数量级
	double magnitudeValue = 1.0;

	if (fNumber == 0.0)
		return(0.0);

	//是否为负数
	bool bNegative;
	bNegative = (fNumber<0) ? true : false;

	double positiveNumber = abs(fNumber);
	if (positiveNumber == 1)
	{//等于1	
		magnitudeValue = 1.0;
	}
	else if (positiveNumber<1.0)
	{//小于1
		while (positiveNumber<1.0)
		{
			positiveNumber *= 10.0;
			magnitudeValue /= 10.0;
		}
	}
	else
	{//大于1
		while (positiveNumber>1.0)
		{
			positiveNumber /= 10.0;
			magnitudeValue *= 10.0;
		}
		magnitudeValue /= 10.0;
	}
	return magnitudeValue;
}
float Jie_Ds::FindStep(float StepMin, bool bUporDown)
{
	string str;
	float dStep, dStepOld, RetVal;
	if (StepMin == 0.0f)
		return(0.0f);

	//是否为负数
	bool bNegative = false;
	bNegative = (StepMin<0) ? true : false;

	if (bNegative)
	{//负数的绝对值的规整方向与之规整方向相反
		bUporDown = !bUporDown;
	}

	//首先按照正数计算
	StepMin = (float)fabs(StepMin);     //若小于0,则取绝对值

	if (!bUporDown)
	{//向下规整
	 //计算渐变母值
		dStep = (float)GetMagnitude(double(StepMin));//得到这个数的数量级
		dStep *= 10;

		while (StepMin < dStep)
		{
			dStepOld = dStep;
			dStep = dStep / 2.5f;
			if (StepMin < dStep)
			{
				dStepOld = dStep;
				dStep = dStep / 2;
			}
			if (StepMin < dStep)
			{
				dStepOld = dStep;
				dStep = dStep / 2.0f;
			}
		}
		RetVal = dStep;
	}
	else
	{//向上规整
	 //计算渐变母值
		dStep = (float)GetMagnitude(double(StepMin));//得到这个数的数量级
		while (StepMin > dStep)
		{
			dStep *= 2.5f;
			if (StepMin > dStep)
			{
				dStep *= 2.0f;
			}
			else
			{
				dStep /= 2.5f;
				dStep *= 2.0f;
			}
		}
		RetVal = dStep;
		//	RetVal = FindMaxFloatStep(StepMin);
	}

	//负数按正数规整后还原其符号
	if (bNegative)
	{
		RetVal *= -1;
	}
	return(RetVal);
}
void Jie_Ds::CalcSameArray()
{
	//等值线最小显示值
	m_Show_MinValue = FindStep((float)m_ZMin, false);//向下规整小值

	//等值线最大显示值
	m_Show_MaxValue = FindStep((float)(m_ZMax), true);//向上规整大值

																//等值线间隔
	m_Show_JianGeValue = (m_Show_MaxValue - m_Show_MinValue) / 2;

	m_sameArray.RemoveAll();
	double minvalue = m_Show_MinValue;
	while (minvalue < m_Show_MaxValue)
	{
		m_sameArray.Add(minvalue);
		minvalue += m_Show_JianGeValue;
	}

	if (minvalue != m_Show_MaxValue)
	{
		m_sameArray.Add(m_Show_MaxValue);
	}
}
void Jie_Ds::SetTrackValue(vector<double> Track)
{
	m_TrackValue.clear();
	for (int i = 0; i < (int)Track.size(); i++)
	{
		m_TrackValue.push_back(Track[i]);
	}
	/*if ((int)m_TrackValue.size() > 2)
	{
		m_valuedis = abs(m_TrackValue[1] - m_TrackValue[0]);
	}
	else
	{
		m_valuedis = 1.0;
	}*/
}
//得到所有的等值点
//void Jie_Ds::EquivalentPoints(double Value, vector<THRVALUE>&Jie_VirtualIJK)
void Jie_Ds::EquivalentPoints(double Value, vector<TWOVALUE>&VirtualIJ)
{
	IsoLine temp_Line;
	VirtualIJ.clear();
	for (int i = 1; i <= m_XNum; i++)//先求纵边上的等值点
	{
		for (int j = 1; j < m_YNum; j++)
		{
			//if (m_YFault[i][j])
			//{
			//	//并且两边有交点 20130913
			//	Y_FalutEquivalent(i, j, Value, VirtualIJ);//得到特殊的等值点 20131015
			//	continue;
			//}
			//判断此边上是否有等值点
			if ((m_GridPoint[i][j].Z - Value) * (m_GridPoint[i][j + 1].Z - Value) < 0)
			{
				TWOVALUE P;
				P.X = i;
				//这里的t值一定是正的！！！！！
				double t = ((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i][j + 1].Z - m_GridPoint[i][j].Z));
				P.Y = j + t;
				//if (t >= 1.0)
				//{
				//	int sk = 1;
				//	ASSERT(sk);
				//	CString str;
				//	str.Format("xt错误 %f,%f,%f", Value, m_GridPoint[i][j].Z, m_GridPoint[i][j + 1].Z);
				//	//AfxMessageBox(str);
				//}
				TWOVALUE RealyPoint;
				RealyPoint.X = m_GridPoint[i][j].X + t*(m_GridPoint[i][j + 1].X - m_GridPoint[i][j].X);
				RealyPoint.Y = m_GridPoint[i][j].Y + t*(m_GridPoint[i][j + 1].Y - m_GridPoint[i][j].Y);
				/*if (IsInside(RealyPoint))
				{
					continue;
				}
				else*/
				VirtualIJ.push_back(P);
				//Jie_VirtualIJK.push_back(P1);
					temp_Line.Logic.push_back(RealyPoint);
					//temp_Line.Value = Value;
					//Jie_IsoLine.push_back(temp_Line);
					//VirtualIJ1.push_back(RealyPoint);
			}
		}
	}
	for (int j = 1; j <= m_YNum; j++)//求横边上的等值点
	{
		for (int i = 1; i < m_XNum; i++)
		{
			//if (m_XFault[i][j])
			//{
			//	//并且两边有交点 20130913
			//	X_FalutEquivalent(i, j, Value, VirtualIJ);//得到特殊的等值点 20131015
			//	continue;
			//}
			if ((m_GridPoint[i][j].Z - Value) * (m_GridPoint[i + 1][j].Z - Value) < 0)
			{
				TWOVALUE P;
				double t = ((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i + 1][j].Z - m_GridPoint[i][j].Z));
				P.X = i + t;
				P.Y = j;
				if (t >= 1.0 || t<0.0)
				{
					//AfxMessageBox("yt错误");
				}
				TWOVALUE RealyPoint;
				RealyPoint.X = m_GridPoint[i][j].X + t*(m_GridPoint[i + 1][j].X - m_GridPoint[i][j].X);
				RealyPoint.Y = m_GridPoint[i][j].Y + t*(m_GridPoint[i + 1][j].Y - m_GridPoint[i][j].Y);
				/*if (IsInside(RealyPoint))
				{
					continue;
				}*/
				VirtualIJ.push_back(P);
				//Jie_VirtualIJK.push_back(P1);
				//VirtualIJ1.push_back(RealyPoint);
				temp_Line.Logic.push_back(RealyPoint);
				
			}
		}
	}
	temp_Line.Value = Value;
	Jie_IsoLine.push_back(temp_Line);
	//std::cout << "================================================" << endl;
	//std::cout << "================================================" << endl;
	//for (int i = 0; i < VirtualIJ1.size(); i++)
	//{
	//	std::cout << "v " << VirtualIJ1[i].X << " " << VirtualIJ1[i].Y << " " << "0" << endl;
    //}
	//std::cout << "================================================" << endl;
	//std::cout << "================================================" << endl;
}

//这里要注意，这种通过vector<vector>定义的二维数组，第一个值在左上角，而非左下角
void Jie_Ds::SignBorder(vector<vector<double>>& biaoji1, vector<vector<double>>& biaoji2)
{
	biaoji1.clear();
	biaoji2.clear();
	biaoji1.resize(m_XNum + 1);
	biaoji2.resize(m_XNum + 2);

	for (int i = 0; i <= m_XNum; i++)
	{
		biaoji1[i].resize(m_YNum + 2);
		biaoji2[i].resize(m_YNum + 1);
		biaoji1[i][0] = 1;
		biaoji1[i][m_YNum + 1] = 1;
		biaoji2[i][0] = 1;
		biaoji2[i][m_YNum] = 1;
	}
	
	biaoji2[m_XNum + 1].resize(m_YNum + 1);
	biaoji2[m_XNum + 1][0] = 1;
	biaoji2[m_XNum + 1][m_YNum] = 1;

	for (int i = 1; i<m_YNum; i++)
	{
		biaoji2[0][i] = 1;
		biaoji2[m_XNum + 1][i] = 1;
		biaoji1[0][i] = 1;
		biaoji1[m_XNum][i] = 1;
	}
	biaoji1[0][m_YNum] = 1;
	biaoji1[m_XNum][m_YNum] = 1;
}
//转化一个点
void Jie_Ds::GetReallyPoint(TWOVALUE A, TWOVALUE &B)
{
	double X1 = A.X;
	double Y1 = A.Y;
	int X0 = (int)X1;
	int Y0 = (int)Y1;
	double Mark = 0, C, D;
	//说明是一个网格点
	if (X1 - X0 == 0 && Y1 - Y0 == 0)
	{
		C = m_GridPoint[X0][Y0].X;
		D = m_GridPoint[X0][Y0].Y;
	}
	//说明在横轴
	else if (X1 - X0 == 0 && Y1 - Y0 != 0)
	{
		Mark = abs(Y1 - Y0);
		C = (1 - Mark) * m_GridPoint[X0][Y0].X + Mark * (m_GridPoint[X0][Y0 + 1].X);//类似于线性插值
		D = (1 - Mark) * m_GridPoint[X0][Y0].Y + Mark * (m_GridPoint[X0][Y0 + 1].Y);
	}
	else if (X1 - X0 != 0 && Y1 - Y0 == 0)
	{
		Mark = abs(X1 - X0);
		C = (1 - Mark) * m_GridPoint[X0][Y0].X + Mark * (m_GridPoint[X0 + 1][Y0].X);
		D = (1 - Mark) * m_GridPoint[X0][Y0].Y + Mark * (m_GridPoint[X0 + 1][Y0].Y);
	}
	B.X = C;
	B.Y = D;

}
void Jie_Ds::TrackRight(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{
	int Num = (int)Line.Logic.size();
	int i = (int)Line.Logic[Num - 1].X;
	int j = (int)Line.Logic[Num - 1].Y;
	double A = m_GridPoint[i][j].Z;
	double B = m_GridPoint[i + 1][j].Z;
	double C = m_GridPoint[i + 1][j + 1].Z;
	double D = m_GridPoint[i][j + 1].Z;
	double Z = (A + B + C + D) / 4;
	//if (m_NoFault[i][j]>-1)					//网格内有断层
	//{
	//	TrackRightFault(Value, Line, IsoReal, X_Used, Y_Used);//追踪网格内
	//	return;
	//}
	//正常追踪
	if ((A - Value) * (B - Value) < 0 && (D - Value) * (C - Value) < 0) //追踪三个等值点
	{
		if (X_Used[i][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - A) / (A - B));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪下边
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j + 1] == 0 && ((D - Value) * (Z - Value) < 0 || ((B - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪上边
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}

	}
	else  //追踪一个等值点
	{
		if ((A - Value) * (B - Value) < 0 && X_Used[i][j] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - A) / (A - B));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪下边
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if ((D - Value) * (C - Value) < 0 && X_Used[i][j + 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪上边
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && (B - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//向右追踪
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackDown(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{
	int Num = (int)Line.Logic.size();
	int i = (int)Line.Logic[Num - 1].X;
	int j = (int)Line.Logic[Num - 1].Y;
	double A = m_GridPoint[i][j].Z;
	double B = m_GridPoint[i + 1][j].Z;
	double C = m_GridPoint[i + 1][j - 1].Z;
	double D = m_GridPoint[i][j - 1].Z;
	double Z = (A + B + C + D) / 4;
	//if (m_NoFault[i][j - 1]>-1)					//网格内有断层
	//{
	//	TrackDownFault(Value, Line, IsoReal, X_Used, Y_Used);//追踪网格内
	//	return;
	//}
	//正常追踪
	if ((A - Value) * (D - Value) < 0 && (B - Value) * (C - Value) < 0) //追踪三个等值点
	{
		if (Y_Used[i][j - 1] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j - abs((Value - A) / (D - A));//这里相当于tmpV.Y = j -1+ abs((Value - D) / (A - D));因为-1+abs((Value - D) / (A - D))=-abs((Value - A) / (D - A))
			Line.Logic.push_back(tmpV);
			Y_Used[i][j - 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪左边
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j - 1] == 0 && ((B - Value) * (Z - Value) < 0 || ((D - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j - abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j - 1] = 1;

			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);
			//追踪右边
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else  //追踪一个等值点
	{
		if ((A - Value) * (D - Value) < 0 && Y_Used[i][j - 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j - abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j - 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪左边
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if ((B - Value) * (C - Value) < 0 && Y_Used[i + 1][j - 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j - abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j - 1] = 1;

			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);
			//追踪右边
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j - 1] == 0 && (D - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j - 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j - 1] = 1;
			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);
			//向下追踪
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackLeft(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{
	int Num = (int)Line.Logic.size();
	int i = (int)Line.Logic[Num - 1].X;
	int j = (int)Line.Logic[Num - 1].Y;
	double A = m_GridPoint[i][j].Z;
	double B = m_GridPoint[i - 1][j].Z;
	double C = m_GridPoint[i - 1][j + 1].Z;
	double D = m_GridPoint[i][j + 1].Z;
	double Z = (A + B + C + D) / 4;
	//if (m_NoFault[i - 1][j]>-1)					//网格内有断层
	//{
	//	TrackLeftFault(Value, Line, IsoReal, X_Used, Y_Used);//追踪网格内
	//	return;
	//}
	//正常追踪
	if ((A - Value) * (B - Value) < 0 && (D - Value) * (C - Value) < 0) //追踪三个点
	{
		if (X_Used[i - 1][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - A) / (B - A));//这里相当于tmpV.X = i -1+ abs((Value - B) / (A - B));因为-1+abs((Value - B) / (A - B))=-abs((Value - A) / (B - A))
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪下边
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i - 1][j + 1] == 0 && ((B - Value) * (Z - Value) < 0 || ((D - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j + 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪上边
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else    //追踪一个点
	{
		if ((A - Value) * (B - Value) < 0 && X_Used[i - 1][j] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - A) / (B - A));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪下边
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if ((D - Value) * (C - Value) < 0 && X_Used[i - 1][j + 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j + 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪上边
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i - 1][j] == 0 && (B - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i - 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//向右追踪
			//TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
			//修改：向左追踪
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackUp(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{
	int Num = (int)Line.Logic.size();
	//当前最后一个等值点的虚拟坐标
	int i = (int)Line.Logic[Num - 1].X;
	int j = (int)Line.Logic[Num - 1].Y;
	double A = m_GridPoint[i][j].Z;
	double B = m_GridPoint[i + 1][j].Z;
	double C = m_GridPoint[i + 1][j + 1].Z;
	double D = m_GridPoint[i][j + 1].Z;
	double Z = (A + B + C + D) / 4;
	//if (m_NoFault[i][j]>-1)					//网格内有断层
	//{
	//	TrackUPFault(Value, Line, IsoReal, X_Used, Y_Used);//追踪网格内
	//	return;
	//}

	if ((A - Value) * (D - Value) < 0 && (B - Value) * (C - Value) < 0)//追踪三个等值点
	{
		if (Y_Used[i][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j + abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪左边
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && ((B - Value) * (Z - Value) < 0 || ((D - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪右边
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else	//只能追踪一个等值点 
	{

		if (Y_Used[i][j] == 0 && (A - Value) * (D - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j + abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪左边
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && (B - Value) * (C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//追踪右边
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j + 1] == 0 && (D - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//实际值
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//向上追踪
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackX(TWOVALUE A, double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{

	int i = (int)A.X;
	int j = (int)A.Y;
	if (j > 1 && j < m_YNum)					//非边界点,处于第一行与最后一行之间
	{
		//虚拟值
		IsoLine LineA;
		IsoLine LineB;
		LineA.Logic.push_back(A);
		LineB.Logic.push_back(A);


		IsoLine RLineA;
		IsoLine RLineB;


		//实际值
		TWOVALUE tmpR;
		//得到真实的等值点的坐标值
		GetReallyPoint(A, tmpR);
		RLineA.Logic.push_back(tmpR);
		RLineB.Logic.push_back(tmpR);
		TrackUp(Value, LineA, RLineA, X_Used, Y_Used);

		int Num = RLineA.Logic.size() - 1;
		if (RLineA.Logic[0].X == RLineA.Logic[Num].X && RLineA.Logic[0].Y == RLineA.Logic[Num].Y)
		{
			IsoReal = RLineA;
			X_Used[i][j] = 1;
			return;
		}
		X_Used[i][j] = 1;
		TrackDown(Value, LineB, RLineB, X_Used, Y_Used);

		//TODO:将两数组联系起来
		int NumA = (int)LineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //不加入首点
			{
				TWOVALUE PointA = LineA.Logic[k];
				Line.Logic.push_back(PointA);
			}
		}
		for (int k = 0; k < (int)LineB.Logic.size(); k++)
		{
			TWOVALUE PointA = LineB.Logic[k];
			Line.Logic.push_back(PointA);
		}

		//实际坐标相连
		NumA = (int)RLineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA-1; k > 0; k--) //不加入首点
			{
				TWOVALUE PointA = RLineA.Logic[k];
				IsoReal.Logic.push_back(PointA);
			}
		}
		for (int k = 0; k < (int)RLineB.Logic.size(); k++)
		{
			TWOVALUE PointA = RLineB.Logic[k];
			IsoReal.Logic.push_back(PointA);
		}

	}
	else if (j == 1)//第一行
	{

		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		X_Used[i][j] = 1;
	}
	else if (j == m_YNum)//最后一行
	{

		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		X_Used[i][j] = 1;
	}
	//X_Used[i][j] = 1;
}
void Jie_Ds::TrackY(TWOVALUE A, double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{

	int i = (int)A.X;
	int j = (int)A.Y;
	if (i > 1 && i < m_XNum)					//非边界点，处于第一列与最后一列之间
	{
		IsoLine LineA;
		IsoLine LineB;
		LineA.Logic.push_back(A);
		LineB.Logic.push_back(A);

		IsoLine RLineA;
		IsoLine RLineB;

		//实际值
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		RLineA.Logic.push_back(tmpR);
		RLineB.Logic.push_back(tmpR);
		TrackLeft(Value, LineA, RLineA, X_Used, Y_Used);
		int Num = RLineA.Logic.size() - 1;
		if (RLineA.Logic[0].X == RLineA.Logic[Num].X && RLineA.Logic[0].Y == RLineA.Logic[Num].Y)
		{
			IsoReal = RLineA;
			Y_Used[i][j] = 1;
			return;
		}
		Y_Used[i][j] = 1;
		TrackRight(Value, LineB, RLineB, X_Used, Y_Used);

		//虚拟坐标相连
		int NumA = (int)LineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //不加入首点
			{
				TWOVALUE PointA = LineA.Logic[k];
				Line.Logic.push_back(PointA);
			}
		}
		for (int k = 0; k < (int)LineB.Logic.size(); k++)
		{
			TWOVALUE PointA = LineB.Logic[k];
			Line.Logic.push_back(PointA);
		}
		//实际坐标相连
		NumA = (int)RLineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //不加入首点
			{
				TWOVALUE PointA = RLineA.Logic[k];
				IsoReal.Logic.push_back(PointA);
			}
		}
		for (int k = 0; k < (int)RLineB.Logic.size(); k++)
		{
			TWOVALUE PointA = RLineB.Logic[k];
			IsoReal.Logic.push_back(PointA);
		}
	}
	else if (i == 1)//第一列
	{
		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		Y_Used[i][j] = 1;
	}
	else if (i == m_XNum)  //最后一列
	{
		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		Y_Used[i][j] = 1;
	}
	//Y_Used[i][j] = 1;
}
//从初始值追踪等值点
void Jie_Ds::TrackPoint(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used, TWOVALUE First)
{
	//这里之所以要进行这一步判断，是因为这里相当于一个取整的操作
	//double t = abs((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i][j + 1].Z - m_GridPoint[i][j].Z));
	//double t = abs((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i+1][j].Z - m_GridPoint[i][j].Z));
	//这里的t相当在网格边上线性插值前面的系数，经过验证t在这里一定是正的
	int i = int(First.X);
	int j = (int)First.Y;
	//说明等值点是网格点
	if (i == First.X && j == First.Y)
	{
		//AfxMessageBox("等值点是网格点");			//用于测试
	}

	if (j == First.Y)								//横轴
	{
		if (X_Used[i][j] == 0)						//没有被使用过
		{
			//cout << "hello1" << endl;
			//TODO:追踪横轴
			//X_Used[int(First.X)][j] = 1;			//首点不作判断，方便以后追踪到首点
			TrackX(First, Value, Line, IsoReal, X_Used, Y_Used);	//追踪横边
		}
	}
	else											//纵轴
	{
		if (Y_Used[i][j] == 0)						//没有被使用过
		{
			//cout << i << " " << j << endl;
			//cout << "hello2" << endl;
			//TODO:追踪纵轴
			//Y_Used[int(First.X)][j] = 1;			//首点不作判断，方便以后追踪到首点
			TrackY(First, Value, Line, IsoReal, X_Used, Y_Used);	//追踪纵边
		}
	}
}
void Jie_Ds::TrackOneValue(double Value)
{
	vector<TWOVALUE> VirtualIJ;
	vector<THRVALUE> Jie_VirtualIJK;
	EquivalentPoints(Value, VirtualIJ);//得到所有的等值点
	vector<vector<double>> X_Used;
	vector<vector<double>> Y_Used;
    //给标记赋值,这里要注意这两个数组的遍历方式
	SignBorder(X_Used, Y_Used);						
	//cout << "======================================================" << endl;
	//cout << "以下是关于标记数组的输出信息：" << endl;
	//for (int i = 0; i < X_Used.size(); i++)
	//{
	//	for (int j = 0; j < X_Used[i].size(); j++)
	//	{
	//		cout << X_Used[i][j] << " ";
	//		//flag_x[i][j] = X_Used[i][j];
	//	}
	//	cout << endl;
	//}
	//cout << endl;
	//for (int i = 0; i < Y_Used.size(); i++)
	//{
	//	for (int j = 0; j < Y_Used[i].size(); j++)
	//	{
	//		cout << Y_Used[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	//cout << "======================================================" << endl;													
	//等值线追踪
	//cout << (int)VirtualIJ.size() << " ";
	for (int i = 0; i < (int)VirtualIJ.size(); i++)
	{
		IsoLine ISOLineOne;	//存储虚拟坐标
		IsoLine IsoReal;	//存储实际坐标
		TWOVALUE First = VirtualIJ[i];
		TrackPoint(Value, ISOLineOne, IsoReal, X_Used, Y_Used, First);
		if (IsoReal.Logic.size() > 20)
		{
			ISOLineOne.Value = Value;
			IsoReal.Value = Value;
			m_IsoLine.push_back(ISOLineOne);
			m_IsoRealLine.push_back(IsoReal); //20131008测试断层内部的网格线
		}
	}
	/*cout << "======================================================" << endl;
	cout << "以下是执行后关于标记数组的输出信息：" << endl;
	for (int i = 0; i < X_Used.size(); i++)
	{
		for (int j = 0; j < X_Used[i].size(); j++)
		{
			cout << X_Used[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	for (int i = 0; i < Y_Used.size(); i++)
	{
		for (int j = 0; j < Y_Used[i].size(); j++)
		{
			cout << Y_Used[i][j] << " ";
		}
		cout << endl;
	}
	cout << "======================================================" << endl;*/
	
}
void Jie_Ds::IsolineTracking()
{

	//在此处清空等值线上数据
	m_IsoLine.clear();
	m_IsoRealLine.clear();	//20131111
	for (int i = 0; i < (int)m_TrackValue.size(); i++)
	{
		double Value = m_TrackValue[i];
		if (Value > m_ZMax || Value < m_ZMin)
		{
			continue;
		}
		else
		{
			TrackOneValue(Value);
		}
		std::ofstream out300("current_isoline.obj");
		for (int i = 0; i < m_IsoRealLine.size(); i++)
		{
			for (int j = 0; j < m_IsoRealLine[i].Logic.size(); j++)
				out300 << "v " << fixed << setprecision(5) << m_IsoRealLine[i].Logic[j].X << " " << fixed << setprecision(5) << m_IsoRealLine[i].Logic[j].Y << " " << "0" << endl;
		}
	}
}
void Jie_Ds::CreateIsoLine()
{
	if (m_ZMax != m_ZMin)
	{
		vector<HwIsoLine> vec;
		m_lsoLines.swap(vec);

		int iso_size = int(m_IsoLine.size());


		for (int i = 0; i < iso_size; i++)
		{
			IsoLine il = m_IsoLine[i];

			HwIsoLine line;
			line.zValue = il.Value;
			line.ptn = il.Logic.size();
			line.pts = new HwIsoPt[line.ptn];

			for (int k = 0; k < il.Logic.size(); k++)
			{
				HwIsoPt hpt;
				hpt.x = il.Logic[k].X;
				hpt.y = il.Logic[k].Y;

				line.pts[k].x = hpt.x;
				line.pts[k].y = hpt.y;
			}

			m_lsoLines.push_back(line);
		}

	}

}



void Jie_Ds::SetOriBoder(vector<TWOVALUE> OriBoder)
{
	Jie_OriBoder.clear();
	for (int i = 0; i < (int)OriBoder.size(); i++)
	{
		Jie_OriBoder.push_back(OriBoder[i]);
	}
}
void Jie_Ds::SetOriISOLine(vector<IsoLine> IsoRealLine)
{
	//for (int i = 0; i < IsoRealLine.size(); i++)
	//	cout << IsoRealLine[i].Logic[0].X << " ";





	//IsoLine L;
	//Jie_IsoRealLine.clear();
	for (int i = 0; i < (int)IsoRealLine.size(); i++)
	{
		/*L.direction0 = IsoRealLine[i].direction0;
		L.direction1 = IsoRealLine[i].direction1;
		L.Index = IsoRealLine[i].Index;
		cout << IsoRealLine[i].Logic.size();
		for (int j = 0; j < IsoRealLine[i].Logic.size(); j++)
		{
			cout << IsoRealLine[i].Logic.size();
			cout << IsoRealLine[i].Logic[j].X << " ";
		}
			
		
		L.Value = IsoRealLine[i].Value;
		cout << "value = ;" << L.Value << endl;*/
		Jie_IsoRealLine.push_back(IsoRealLine[i]);
	}
}
bool Jie_Ds::IsinLineK(TWOVALUE Star, TWOVALUE End, TWOVALUE Point)
{
	End.X = End.X - Star.X;
	End.Y = End.Y - Star.Y;
	Point.X = Point.X - Star.X;
	Point.Y = Point.Y - Star.Y;
	Star.X = 0;
	Star.Y = 0;					//平移，对判断是否一条直线不造成影响··将大数化为小数~能提高计算精度
	if (Point.X *(Point.X - End.X) >0 && Point.Y *(Point.Y - End.Y) > 0)
	{
		return false;
	}
	double A = Star.Y - End.Y;
	double B = End.X - Star.X;
	//double C = -A * Star.X- B * Star.Y;
	double C = 0;
	double D = sqrt(pow(A, 2) + pow(B, 2));
	double d = A * (Point.X - Star.X) + B * (Point.Y - Star.Y);
	d = d / D;
	if (abs(d)<0.00001)
	{
		return true;
	}
	return false;

}
//线段是否相交,index虚拟序号
bool Jie_Ds::ISIntersect(TWOVALUE p1_s, TWOVALUE p1_e, TWOVALUE p2_s, TWOVALUE p2_e, double &index)
{
	double Xmax_1, Xmax_2, Xmin_1, Xmin_2, Ymax_1, Ymax_2, Ymin_1, Ymin_2;
	double V1, V2, V3, V4;

	if (p1_s.X > p1_e.X)
	{
		Xmax_1 = p1_s.X;
		Xmin_1 = p1_e.X;
	}
	else
	{
		Xmax_1 = p1_e.X;
		Xmin_1 = p1_s.X;
	}


	if (p1_s.Y > p1_e.Y)
	{
		Ymax_1 = p1_s.Y;
		Ymin_1 = p1_e.Y;
	}
	else
	{
		Ymax_1 = p1_e.Y;
		Ymin_1 = p1_s.Y;
	}


	if (p2_s.X > p2_e.X)
	{
		Xmax_2 = p2_s.X;
		Xmin_2 = p2_e.X;
	}
	else
	{
		Xmax_2 = p2_e.X;
		Xmin_2 = p2_s.X;
	}

	if (p2_s.Y > p2_e.Y)
	{
		Ymax_2 = p2_s.Y;
		Ymin_2 = p2_e.Y;
	}
	else
	{
		Ymax_2 = p2_e.Y;
		Ymin_2 = p2_s.Y;
	}


	if (Xmax_1 < Xmin_2 || Xmin_1 > Xmax_2 || Ymin_1 > Ymax_2 || Ymax_1 < Ymin_2)   //两线段最小矩形不相交，得出两线段不相交
		return FALSE;
	else
	{
		V1 = (p1_e.X - p1_s.X) * (p2_s.Y - p1_s.Y) - (p2_s.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V2 = (p1_e.X - p1_s.X) * (p2_e.Y - p1_s.Y) - (p2_e.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V3 = (p2_e.X - p2_s.X) * (p1_s.Y - p2_s.Y) - (p1_s.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		V4 = (p2_e.X - p2_s.X) * (p1_e.Y - p2_s.Y) - (p1_e.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		//特殊情况点在直线上时候   20131015
		BOOL Mark = FALSE;
		TWOVALUE TempPoint;
		if (abs(V1)< 0.00000001)
		{
			TempPoint = p2_s;
			Mark = TRUE;
		}
		else if (abs(V2)< 0.00000001)
		{
			TempPoint = p2_e;
			Mark = TRUE;
		}
		else if (abs(V3)< 0.00000001)
		{
			TempPoint = p1_s;
			Mark = TRUE;
		}
		else if (abs(V4)< 0.00000001)
		{
			TempPoint = p1_e;
			Mark = TRUE;
		}
		if (Mark)
		{
			double t = 0;
			if (abs(TempPoint.X - p1_s.X)<abs(TempPoint.Y - p1_s.Y))
			{
				t = (TempPoint.Y - p1_s.Y) / (p1_e.Y - p1_s.Y);/////////////////////
			}
			else
			{
				t = (TempPoint.X - p1_s.X) / (p1_e.X - p1_s.X);
			}
			index = index + t;
			return TRUE;
		}
		if (V3 * V4 <= 0 && V1 * V2 <= 0)
		{
			double dy = p1_e.Y - p1_s.Y;
			double dx = p1_e.X - p1_s.X;
			double t = 0;
			if (abs(p2_e.X - p2_s.X)<0.00000000001)
			{
				t = (p2_e.X - p1_s.X) / dx;/////////////////////
			}
			else
			{
				double k = (p2_s.Y - p2_e.Y) / (p2_s.X - p2_e.X);
				t = ((p2_s.Y - p1_s.Y) - k*p2_s.X + k*p1_s.X) / (dy - k*dx);
			}
			index = index + t;
			return TRUE;
		}
		else
			return FALSE;
	}
}
//点A是否在Line内部
bool Jie_Ds::IsInside(TWOVALUE A, vector<TWOVALUE> Line)
{
	//cout << "hhh" << endl;
	for (int i = 0; i < (int)Line.size() - 1; i++)
	{
		if (IsinLineK(Line[i], Line[i + 1], A))		//点在边界上  20131023
		{
			return TRUE;
		}
	}
	double Xmax = Line[0].X;
	double Xmin = Line[0].X;
	double Ymax = Line[0].Y;
	double Ymin = Line[0].Y;
	for (int i = 1; i < (int)Line.size(); i++)
	{
		if (Line[i].X > Xmax)
		{
			Xmax = Line[i].X;
		}
		else if (Line[i].X < Xmin)
		{
			Xmin = Line[i].X;
		}

		if (Line[i].Y > Ymax)
		{
			Ymax = Line[i].Y;
		}
		else if (Line[i].Y < Ymin)
		{
			Ymin = Line[i].Y;
		}
	}
	if (A.X > Xmax || A.X <Xmin || A.Y>Ymax || A.Y < Ymin)
	{
		return FALSE;
	}
	//射线判断是否在内部
	TWOVALUE B = A;
	B.X = m_XMin - 200.0;     //AB射线
	int Sum = 0;
	for (int i = 0; i < (int)Line.size() - 1; i++)
	{
		double t = 0;
		TWOVALUE LineA = Line[i];
		TWOVALUE LineB = Line[i + 1];
		if (ISIntersect(LineA, LineB, A, B, t))
		{
			Sum = Sum + 1;
		}
	}
	if (Sum % 2 == 1)//交点为奇数个则在内部
	{
		return TRUE;
	}
	return FALSE;
}
bool Jie_Ds::L2L_Intersect(TWOVALUE p1_s, TWOVALUE p1_e, TWOVALUE p2_s, TWOVALUE p2_e, TWOVALUE &Point)
{
	double Xmax_1, Xmax_2, Xmin_1, Xmin_2, Ymax_1, Ymax_2, Ymin_1, Ymin_2;
	double V1, V2, V3, V4;

	if (p1_s.X > p1_e.X)
	{
		Xmax_1 = p1_s.X;
		Xmin_1 = p1_e.X;
	}
	else
	{
		Xmax_1 = p1_e.X;
		Xmin_1 = p1_s.X;
	}


	if (p1_s.Y > p1_e.Y)
	{
		Ymax_1 = p1_s.Y;
		Ymin_1 = p1_e.Y;
	}
	else
	{
		Ymax_1 = p1_e.Y;
		Ymin_1 = p1_s.Y;
	}


	if (p2_s.X > p2_e.X)
	{
		Xmax_2 = p2_s.X;
		Xmin_2 = p2_e.X;
	}
	else
	{
		Xmax_2 = p2_e.X;
		Xmin_2 = p2_s.X;
	}

	if (p2_s.Y > p2_e.Y)
	{
		Ymax_2 = p2_s.Y;
		Ymin_2 = p2_e.Y;
	}
	else
	{
		Ymax_2 = p2_e.Y;
		Ymin_2 = p2_s.Y;
	}


	if (Xmax_1 < Xmin_2 || Xmin_1 > Xmax_2 || Ymin_1 > Ymax_2 || Ymax_1 < Ymin_2)   //两线段最小矩形不相交，得出两线段不相交
		return false;
	else   //利用向量的叉积性质，当其中一条线段的两个端点在另一条线段的同一侧时，不相交。否则，相交。
		   //( Q1 - P1 )×( P2 - P1) * ( P2 - P1)×(Q2 - P1) >= 0。
		   //叉积的计算公式为:  P1 X P2 = x1y2 - x2y1;
	{
		V1 = (p1_e.X - p1_s.X) * (p2_s.Y - p1_s.Y) - (p2_s.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V2 = (p1_e.X - p1_s.X) * (p2_e.Y - p1_s.Y) - (p2_e.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V3 = (p2_e.X - p2_s.X) * (p1_s.Y - p2_s.Y) - (p1_s.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		V4 = (p2_e.X - p2_s.X) * (p1_e.Y - p2_s.Y) - (p1_e.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		//特殊情况点在直线上时候   20131015
		if (abs(V1)< 0.00000001)
		{
			Point = p2_s;
			return true;
		}
		else if (abs(V2)< 0.00000001)
		{
			Point = p2_e;
			return true;
		}
		else if (abs(V3)< 0.00000001)
		{
			Point = p1_s;
			return true;
		}
		else if (abs(V4)< 0.00000001)
		{
			Point = p1_e;
			return true;
		}

		else if (V3 * V4 <= 0 && V1 * V2 <= 0)
		{
			double dy = p1_e.Y - p1_s.Y;
			double dx = p1_e.X - p1_s.X;
			double t = 0;
			if (abs(p2_e.X - p2_s.X)<0.00000000001)
			{
				t = (p2_e.X - p1_s.X) / dx;/////////////////////
				Point.X = p2_s.X;
				Point.Y = t*dy + p1_s.Y;
			}
			else
			{
				double k = (p2_s.Y - p2_e.Y) / (p2_s.X - p2_e.X);
				t = ((p2_s.Y - p1_s.Y) - k*p2_s.X + k*p1_s.X) / (dy - k*dx);
				Point.X = t * dx + p1_s.X;
				Point.Y = t * dy + p1_s.Y;
			}
			return true;
		}
	}
	return false;

}
bool Jie_Ds::GetCrossPt(TWOVALUE Star, TWOVALUE End, vector<TWOVALUE> Line, TWOVALUE &A)
{
	for (int j = 0; j < (int)Jie_OriBoder.size() - 1; j++)
	{
		TWOVALUE C = Jie_OriBoder[j];
		TWOVALUE B = Jie_OriBoder[j + 1];
		if (L2L_Intersect(C, B, Star, End, A))
		{
			return TRUE;
		}
	}
	return FALSE;
}
void Jie_Ds::DleaIso(IsoLine &OneIso, vector<IsoLine> &NewIso)
{
	IsoLine TempLine;
	TempLine.Value = OneIso.Value;
	//cout << "Jie_OriBorder size is :" << Jie_OriBoder.size() << endl;
	//cout << "当前的坐标：" << OneIso.Logic[0].X << " " << OneIso.Logic[0].X << endl;
	if (OneIso.Logic.size() > 0)
	{
		bool flag = IsInside(OneIso.Logic[0], Jie_OriBoder);
		if (!flag)//首点不在原始边界内，删除以前的点，并加入边界点
		{
			for (int i = 1; i < (int)OneIso.Logic.size(); i++)
			{
				if (IsInside(OneIso.Logic[i], Jie_OriBoder))
				{
					TWOVALUE Star = OneIso.Logic[i - 1];
					TWOVALUE End = OneIso.Logic[i];
					TWOVALUE A;
					if (GetCrossPt(Star, End, Jie_OriBoder, A))//得到边界和这两个点构成的线段的交点
					{
						TempLine.Logic.push_back(A);//此交点作为等值线的其中的一点
					}
					else
					{
						int sk = 0;
						//AfxMessageBox("111");
					}
					for (int j = i; j < (int)OneIso.Logic.size(); j++)
					{
						TempLine.Logic.push_back(OneIso.Logic[j]);//插入坐标
					}
					break;
				}
			}
		}
		else
		{
			TempLine = OneIso;
		}
		if ((int)TempLine.Logic.size() <= 1)	//点全部在原始边界外面
		{
			return;
		}
		IsoLine StarLine;
		StarLine.Value = TempLine.Value;
		StarLine.Logic.push_back(TempLine.Logic[0]);
		OneIso.Logic.clear();
		for (int i = 1; i < (int)TempLine.Logic.size(); i++)  //从第二个点开始插入坐标
		{
			if (IsInside(TempLine.Logic[i], Jie_OriBoder))
			{
				StarLine.Logic.push_back(TempLine.Logic[i]);
			}
			else
			{
				//CString str;
				//str.Format("首个边界点%d,%f     %f",i,TempLine.Logic[0].X,TempLine.Logic[0].Y);
				//AfxMessageBox(str);

				TWOVALUE Star = TempLine.Logic[i - 1];
				TWOVALUE End = TempLine.Logic[i];
				TWOVALUE A;
				if (GetCrossPt(Star, End, Jie_OriBoder, A))
				{
					StarLine.Logic.push_back(A);
				}
				else
				{
					int sk = 0;
					//AfxMessageBox("111");
				}

				for (int j = i; j < (int)TempLine.Logic.size(); j++)
				{
					OneIso.Logic.push_back(TempLine.Logic[j]);
				}
				break;
			}
		}
		NewIso.push_back(StarLine);
		if (OneIso.Logic.size() > 0)//如果里面还有数据，说明等值线中间部分也被边界给切分了
		{
			DleaIso(OneIso, NewIso);
		}
	}
}

void Jie_Ds::DleaIso()
{

	vector<IsoLine> NewIso;
	for (int i = 0; i < (int)Jie_IsoRealLine.size(); i++)
	{
		vector<IsoLine> TempLines;
		IsoLine LineOne = Jie_IsoRealLine[i];
	
		DleaIso(LineOne, TempLines);
		for (int j = 0; j < (int)TempLines.size(); j++)
		{
			NewIso.push_back(TempLines[j]);
		}
	}

	Jie_IsoRealLine.clear();
	Jie_IsoRealLine = NewIso;//被边界切割后的等值线
}
vector<double> Jie_Ds::Randnum(int k)
{
	vector<double> temp;
	srand((int)time(0));
	for (int i = 1; i <= k / 2; i++)
	{
		double temp_ = random(10);
		temp.push_back(temp_);
	}
	return temp;
}