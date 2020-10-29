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
	//ѡ����С��x,����y��Ϊ��ʼ��,�洢��m_point[0]��
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
	//�������е�͸õ���ɵ������ͳ�ʼ������pt = ptStartA-ptEndA����ʹ�ø�������xֵΪ��ֵ��֮��ĵļнǣ�ȡһ�����ļн�(��С��cos x = (a*b)/(|a|*|b|))
	//����ĳ�ʼ����ָ��Ӧ����ȡһ��y����������������ɣ�����˵��0,1��
	Point.push_back(Temppoint);
	m_point[StarIndex].mark = 1;	//��ʹ�ù�
	if ((int)m_point.size() == 1)
	{
		return Point;		//��������λ��ֵ 20131108
	}
	vector<double> CosS, NiCosS, Dis, NiDis;
	vector<int> IndexS, NiIndexS;
	//�ҵ��ڶ����㣬����ŵ�m_point[1];
	for (int i = 0; i< Num; i++)
	{
		if (m_point[i].mark == 1)
		{
			continue;	//��ʹ�ù�
		}
		double x1 = 0.0;
		double y1 = 1.0;
		double x2 = m_point[i].X - Point[0].X;
		double y2 = m_point[i].Y - Point[0].Y;
		double d = -x2;		//��˽����������Ҫ���ݷ������ж�(x2,y2)�����ɵ�������(x1,y1)������������˳ʱ�뻹����ʱ�뷽��
							//��С��0������˳ʱ�뷽��������0��������ʱ�뷽��������0������

		double d1 = sqrt(pow(x1, 2) + pow(y1, 2));
		double d2 = sqrt(pow(x2, 2) + pow(y2, 2));

		if (d1 == 0 || d2 == 0)
		{
			continue;
		}

		double cosa = (x1*x2 + y1*y2) / (d1*d2);//����Ӧ���Ǵ���������ֱ���ϵ������нǵ�����ֵ
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
			Dis.push_back(d2);	//�洢����
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
	m_point[StarIndex].mark = 0;	//ɾ���׵㱻ʹ�úۼ�
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
			double d = x1 * y2 - y1 * x2;		//��˽��
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
	double XMin = m_oriData[0].X;//m_oriData�д洢ԭʼ���ݣ���λ��Ϣ
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
	//�Ƿ�Ҫ�ж�m_boder�����ݣ�20131112  ֮ǰ�Ѿ�����ôһ���ĸ�ֵ������m_OriBoder = m_Border;
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
	if (m_Border.size() < 3) // û�����ñ߽磬����һ���߽�
	{
		Data dt;
		dt.X = m_XMin; // ����
		dt.Y = m_YMax;
		m_Border.push_back(dt);

		dt.X = m_XMax; // ����
		dt.Y = m_YMax;
		m_Border.push_back(dt);

		dt.X = m_XMax; // ����
		dt.Y = m_YMin;
		m_Border.push_back(dt);

		dt.X = m_XMin; // ����
		dt.Y = m_YMin;
		m_Border.push_back(dt);
	}


	//m_Border.clear();//20131111
	//m_Border.resize(0);

	//for (int i = 0; i< m_Border.size(); i++)
	//{
	//	m_Border.push_back(m_Border[i]); // ����߽�
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

	for (double t = 0.0; t <= 1.0; t += 0.001 / n)//���������
	{
		for (int i = 1; i < n; i++)
		{
			for (int j = 0; j < n - i; j++)
			{
				if (i == 1)//i=1����֪���ƶ������
				{
					control_point[j].x = (1 - t) * input_vertice[j].x + t * input_vertice[j + 1].x;
					control_point[j].y = (1 - t) * input_vertice[j].y + t * input_vertice[j + 1].y;
					continue;
				}
				else//i!=1����һ�ε����������
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
	//׼����B����
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

					if (r == 1)//���Ƶ�
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
			b_spline.push_back(control_point[j]);//���Ƶ����һ��ĵ㣬��Ϊ��õĵ�
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
	//�����е�       
	for (int i = 0; i < count; i++) {
		int nexti = (i + 1) % count;
		Vector2D p;
		p.x = (originPoint[i].x + originPoint[nexti].x) / 2.0;
		p.y = (originPoint[i].y + originPoint[nexti].y) / 2.0;
		midpoints.push_back(p);
	}

	//ƽ���е�  
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
		//�� originPoint[i]��������   
		double addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;
		double addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;
		extrapoints[extraindex].x = originPoint[i].x + addx;
		extrapoints[extraindex].y = originPoint[i].y + addy;

		int extranexti = (extraindex + 1) % (2 * count);
		extrapoints[extranexti].x = midpoints[i].x + offsetx;
		extrapoints[extranexti].y = midpoints[i].y + offsety;
		//�� originPoint[i]��������   
		addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;
		addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;
		extrapoints[extranexti].x = originPoint[i].x + addx;
		extrapoints[extranexti].y = originPoint[i].y + addy;

	}

	Vector2D controlPoint[4];
	//����4���Ƶ㣬��������������  
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
			//u�Ĳ����������ߵ�����  
			u -= 0.01;
			Vector2D tempP = Vector2D(px, py);
			//�������ߵ�   
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
	//	//�����ı���
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
	//	convexBag.push_back(p4);//���Ͻ�Ϊ��1,1���㣬ʵ���ϵģ�0,0����
	//	convexBag.push_back(p1);
	//}
	//else if ((int)convexBag.size() == 2)
	//{
	//	//�����Ӧ�ľ���
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
	//	convexBag.push_back(p4);		 //���Ͻ�Ϊ��1,1���㣬ʵ���ϵģ�0,0����
	//	convexBag.push_back(p1);        //����Ϊɶ���Ͻ�������

	//}
	//else if ((int)convexBag.size() == 4)
	//{
	//	//���һ����
	//	AddPt(convexBag);
	//}

	//double dx = m_XMax - m_XMin;
	//double dy = m_YMax - m_YMin;
	//if (dy - 8 * dx > 0)
	//{
	//	//�޸�ͼ��Ϊ���Σ�
	//}
	//else if (dx - 8 * dy > 0)
	//{
	//	//�޸�ͼ��Ϊ���Σ�
	//}
}
double Jie_Ds::Angle(Data &p0, const Data &p1, const Data &p2)
{
	//���ڶ���������x�᷽������ƽ��һ����λ��ֵ����һ����
	p0.X = p1.X + 1;
	p0.Y = p1.Y;


	/*std::ofstream out("current_triangle_point.obj");
	out << "v " << p0.X << " " << p0.Y << " " << "0" << endl;
	out << "v " << p1.X << " " << p1.Y << " " << "0" << endl;
	out << "v " << p2.X << " " << p2.Y << " " << "0" << endl;

	out.close();*/

	double A, B, X, cross, angle;
	//���� ����
	A = sqrt((p1.X - p0.X) * (p1.X - p0.X) + (p1.Y - p0.Y) *(p1.Y - p0.Y));
	B = sqrt((p1.X - p2.X) * (p1.X - p2.X) + (p1.Y - p2.Y) *(p1.Y - p2.Y));

	if (A == 0 || B == 0)
		return 0;

	//����ֵ
	X = ((p0.X - p1.X) * (p2.X - p1.X) + (p0.Y - p1.Y) * (p2.Y - p1.Y)) / (A*B);

	//����ж���ʱ���Ƿ����180��
	//a = (x1, y1) b = (x2, y2)
	//a��b = x1y2 - x2y1,�����Ϊ����������b��a����ʱ�뷽��,����b��a��˳ʱ�뷽��,�����Ϊ0����a��b����
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
		if (cross < 0)//˳ʱ��
		{
			//cout << "˳ʱ��" << endl;
			//angle = atan2((p2.Y - p1.Y), (p2.X - p1.X));
			temp = (atan(-X / sqrt(-X * X + 1)) + 2 * atan(1.0));
			angle = p*2  - temp;
		}
		else//��ʱ��
		{
			//cout << "��ʱ��" << endl;
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
	//�������������ʲô����
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
	//����ֵ
	side[0] = FourPart[3];
	side[1] = 0;
	side[2] = FourPart[1];
	side[3] = FourPart[2];

	int count_s0 = 0, count_s1 = 0, count_s2 = 0, count_s3 = 0;
	int j = 1, i;
	while (JudgeMax)
	{
		i = 0;

		//���±߼�����С�����ֵ
		if (side[0] <= FourPart[0])
		{
			count_s0++;
			//���㸳ֵ
			downline[j].X = convexBag[side[0]].X;
			downline[j].Y = convexBag[side[0]].Y;


			//��Ϊ�׵�ʱ����ľ��� Ĭ��Ϊ��(0,0)�ľ��� ��
			//Ϊ�׵�ʱ �����븳ֵ 0
			if (j == 1)
				downDis[j].dis = 0;
			else
			//����˵��˱��׵�ľ���
			downDis[j].dis = sqrt((downline[j].X - downline[j - 1].X) * (downline[j].X - downline[j - 1].X)
				+ (downline[j].Y - downline[j - 1].Y) * (downline[j].Y - downline[j - 1].Y)) + downDis[j - 1].dis;

			

			//���������
			side[0] = side[0] + 1;
			//��¼�жϵ�������
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

		//��߸�ֵ  ͬ��
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
		//�ϱ߸�ֵ ͬ��
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

		//�ұ߸�ֵ ͬ��
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

		//�ߵļ�¼�� ����
		j = j + 1;
		//���û������ i=0 ��˵���Ѿ��������� ���е��ϵĵ� ����ѭ��
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
//����upLine��downline��leftLine��rightLine��Dis����ռ�ܵĳ��ȵİ�ֱ�
//����upDis��downdis��leftDis��rightDis��angle��������ǰһ���㹹�ɵ��߶���x��������ļн�
void Jie_Ds::BordersChar(vector<Data>& upLine, vector<Data>&downline, vector<Data>&leftLine, vector<Data>& rightLine,
	vector<Dis>& upDis, vector<Dis>&downdis, vector<Dis>&leftDis, vector<Dis>&rightDis)
{
	bool JudgeMax = true;
	Data triPoint[3];
	int i, j = 2;
	std::ofstream out22("2.obj");
	while (JudgeMax)
	{
		//һ���жϴ�ѭ���Ƿ��м���� ������ i>0˵����ѭ�������㣬����i=0 ����ѭ��
		//����ֵ
		i = 0;
		//����û�г��� �±� ����м���
		if (j < (int)downline.size())
		{
			//����˵�ռ�˱��ܳ��İٷ���
			//����downline.size()��downdis.size()��1��downsize�±��1��ʼ��downsize[1]=0;
			downdis[j].per = downdis[j].dis / downdis[downline.size() - 1].dis;
			//��ֵ��¼�˵���ǰһ����߶�
			triPoint[1] = downline[j - 1];
			triPoint[2] = downline[j];

			//������߶���x��������ĽǶ�
			downdis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//cout << downdis[j].angle << endl;
			//���� >0˵��������
			i = i + 1;
			
			out22 << "v " << downline[downline.size() - 1].X << " " << downline[downline.size() - 1].Y << " " << "0" << endl;
		}
		
		if (j < (int)leftLine.size())
		{
			//����˵�ռ�˱��ܳ��İٷ���
			leftDis[j].per = leftDis[j].dis / leftDis[leftLine.size() - 1].dis;
			//��ֵ��¼�˵���ǰһ����߶�
			triPoint[1] = leftLine[j - 1];
			triPoint[2] = leftLine[j];

			//������߶���x��������ĽǶ�
			leftDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//���� >0˵��������
			i = i + 1;
			
		}

		if (j < (int)upLine.size())
		{
			//����˵�ռ�˱��ܳ��İٷ���
			upDis[j].per = upDis[j].dis / upDis[upLine.size() - 1].dis;
			//��ֵ��¼�˵���ǰһ����߶�
			triPoint[1] = upLine[j - 1];
			triPoint[2] = upLine[j];

			//������߶���x��������ĽǶ�
			upDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//���� >0˵��������
			i = i + 1;
			
		}

		if (j < (int)rightLine.size())
		{
			//����˵�ռ�˱��ܳ��İٷ���
			rightDis[j].per = rightDis[j].dis / rightDis[rightLine.size() - 1].dis;
			//��ֵ��¼�˵���ǰһ����߶�
			triPoint[1] = rightLine[j - 1];
			triPoint[2] = rightLine[j];

			//������߶���x��������ĽǶ�
			rightDis[j].angle = Angle(triPoint[0], triPoint[1], triPoint[2]);
			//���� >0˵��������
			i = i + 1;
			
		}

		//�߽����������
		j = j + 1;
		//i=0��˵���˴�ѭ��û������ ������ѭ��
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

	//*****************�ָ��㷨
	//ȷ���ĵ�����
	//�洢ӳ���ԭ������
	sim1[1].X = convexBag[FourPart[0]].X;//����Ӧ�ô洢����͹���߽������һ�����λ��

	//��Ӧ��һ���ߵ�x�����ֵ
	sim1[2].X = convexBag[FourPart[3]].X - convexBag[FourPart[0]].X;

	//��Ӧ�ڶ����ߵ�x�����ֵ
	sim1[3].X = convexBag[FourPart[1]].X - convexBag[FourPart[0]].X;

	sim1[4].X = convexBag[FourPart[0]].X - convexBag[FourPart[1]].X + convexBag[FourPart[2]].X - convexBag[FourPart[3]].X;

	//����ͬ�� ��Ӧy����ֵ
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
		//�˷ָ�������ľ���ռ�˱߳��ȵİٷ���
		t = double(j - 1) / (YM - 1);//д��Ӧ��������ɣ���������

									 //����
									 //��˵㵽�˱����ľ���
		sim3[j].X = sim1[3].X * t;
		sim3[j].Y = sim1[3].Y * t;

		//�����˵����������(�洢 �㵽�ߵľ���,���˾���Ϊ�˱ߵİٷ���, _�ʹ˵���ǰһ����߶κ�x��������ĽǶ�)
		for (m = 2; m<(int)leftDis.size(); m++)
		{
			//���ָ�� �ٷ����� ͹������İٷ�����ͬʱ
			if (t == leftDis[m - 1].per) //��͹����Ϊ�ָ��ʱ
			{
				//�Ѵ�͹���㸳ֵ����߷ָ������
				DL[j].X = leftLine[m - 1].X;
				DL[j].Y = leftLine[m - 1].Y;
				//���ҵ��˷ָ���������
				break;
			}

			//���ָ��������˱ߵ�β��ʱ
			else if (t == 1)
			{
				//�Ѵ˱ߵ����һ���㸳ֵ�� �ָ������
				DL[j].X = convexBag[FourPart[1]].X;
				DL[j].Y = convexBag[FourPart[1]].Y;
				//���ҵ�������ѭ��
				break;
			}

			//���˷ָ����׵�ľ���ռ�İٷ��� Ϊ����͹������ʱ ��  t��ĳ���߶����˵�İٷ���֮��
			else if (leftDis[m - 1].per < t && t < leftDis[m].per)
			{
				//�˷ָ������͹�����ϵ�λ�õ���͹����һ�˵����ռ��͹���߳��ȵİٷ���
				Tran = leftDis[leftDis.size() - 1].dis* t - leftDis[m - 1].dis;
				//���ݴ˱ߵĽǶ���˵㵽һ�ε�ľ��� ��ֱ�߲������� �ɼ���˵��(x,y)����
				DL[j].X = leftLine[m - 1].X + Tran * cos(leftDis[m].angle);
				DL[j].Y = leftLine[m - 1].Y + Tran * sin(leftDis[m].angle);
				break;
			}
		}

		if ((int)rightDis.size()>1)
		{
			//ͬ�� �����ұ߽�ķָ�� ��������ֵ
			for (m = 2; m< (int)rightDis.size(); m++)
			{
				if ((1 - t) == rightDis[m - 1].per)  //��Ϊĳ���ڵ�ʱ
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
		//����������=3ʱ, ֻ��һ����
		else
		{
			//�˵�Ϊ��4����
			DR[j].X = convexBag[FourPart[3]].X;
			DR[j].Y = convexBag[FourPart[3]].Y;
		}
	}

	// ���ұ߽��Ϸָ��������

	//�������±߽�ķָ��
	for (i = 1; i <= XM; i++)
	{
		//����˵�ռ�ܱ߳��İٷ���
		s = double(i - 1) / (XM - 1);
		sim2[i].X = sim1[2].X * s;
		sim2[i].Y = sim1[2].Y * s;
		sim4[i].X = sim1[4].X * s;
		sim4[i].Y = sim1[4].Y * s;

		//���ϱ߽�ָ�㸳ֵ
		for (m = 2; m< (int)upDis.size(); m++)
		{
			//���˷ָ��İٷ�����͹����İٷ�����ͬʱ ͹����Ϊ�ָ��
			if (s == upDis[m - 1].per) //��Ϊĳ���ڵ�ʱ
			{
				//���ָ�㸳ֵ
				DU[i].X = upLine[m - 1].X;
				DU[i].Y = upLine[m - 1].Y;
				//����ʱ����ѭ��
				break;
			}

			//Ϊβ��ʱ
			else if (s == 1)
			{
				DU[i].X = convexBag[FourPart[2]].X;
				DU[i].Y = convexBag[FourPart[2]].Y;
				//����ѭ��
				break;

			}

			//��Ϊĳ��͹����ʱ
			else if (upDis[m - 1].per < s && s < upDis[m].per)
			{
				Tran = upDis[upDis.size() - 1].dis * s - upDis[m - 1].dis;
				DU[i].X = upLine[m - 1].X + Tran * cos(upDis[m].angle);
				DU[i].Y = upLine[m - 1].Y + Tran * sin(upDis[m].angle);
				break;
			}
		}

		//����ͬ�� �ɼ����±߽�ָ�������ֵ
		for (m = 2; m < (int)downdis.size(); m++)
		{
			if ((1 - s) == downdis[m - 1].per)  //��Ϊĳ���ڵ�ʱ
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
	//��m�ұ�����һ����λ�󣬹���һ��m���������mm
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
	double **mm = new double *[n];		//���ݱȽϴ������Ƚ϶࣬vectorЧ�ʱ�ָ��ͺܶ�20150325
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
	//ͨ�������б任(����˹��ȥ��)ʹԭ�����Ϊ��λ�����ұߵĵ�λ����ԭ���������
	for (k = 0; k < n - 1; k++)
	{
		/*----------------------------------------*/
		//�ӵ� k �С��� k �п�ʼ�����½�������ѡȡ����ֵ����Ԫ�أ�����ס��Ԫ���ڵ��кź��кţ�
		//��ͨ���н������н���������������Ԫ��λ����.��һ����Ϊȫѡ��Ԫ  20140925
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
			mm[i][k] = 0.0;		//��ֹ�������20140929
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

	//����任����ұߵľ���
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			M[i][j] = mm[i][j + n];
	}

	//������ȫѡ��Ԫ����������¼���С��н�������Ϣ���лָ����ָ���ԭ�����£���ȫѡ��Ԫ�����У�
	//�Ƚ������У��У�����лָ���ԭ�����У��У��������У��У��������ָ���
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
	double b0 = 1;			//�ȼ���б��Ϊһ��ֵ
	double b1 = 0;			//��ؾ�Ϊһ��ֵ
	double RR = m_B;		//������ͺ��
	int n = (int)suV.size();
	vector<vector<double>> Arr(n + 1);

	//********������ֵ
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
				d = RR;			//20140422����ͺ����
			}
			Arr[i][j] = abs(b0 * (d)-RR);
		}
	}

	Arr[0][n] = 0;
	Arr[0][0] = 1;

	/************************************************/
	//CString Str = "zֵ \n";
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
	//	str.Format("%f,    %f,    %f��   %f     ��%f\n",suV[i].X,suV[i].Y,suV[i].Z,m_ZMax,m_ZMin);
	//	Str += str;
	//}
	//AfxMessageBox(Str);
	/************************************************/
	//�������������
	Inv(Arr);//�������������Ҫ��ʹ��eigen��Ļ���ֻ��Ҫһ�д���

			 //��ⷽ�̸�
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
	double ZZ = m_B;//ֱ���ͱ��캯��Bֵ
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
			distanceVal = ZZ;	//20140422����ͺ����
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
		//d = pow(d,1.5);		//�Ӵ���뷴�ȷ��Ĵ������Ա���ý��ĵ�Ӱ��������20131217
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
		//	d = d * 1000;								//���Ŷϲ�����СȨֵ20140904��
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
	if ((int)m_oriData.size() < 800)//20150323,800�ھ�����ʹ�ÿ���𷽷����ٶ�4s���ҿ��Խ���
	{
		m_IsK = true;				//��������ʹ�ÿ�����ֵ
	}
	else
	{
		m_IsK = false;				//����̫��ʹ���ݾ��뷴�ȹ�ֵ
	}

	if ((int)m_oriData.size() <= 1)
	{
		return;		//����ֻ��һ����λ������20131107
	}
	if (m_IsK)			//�Ƿ�ʹ�ÿ�����㷨��ֵ��׼ȷ�����ٶ�����201308163
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
				//���ƴ˵������ֵ    û�жϲ�Ĺ���ֵ + ����ֵ
				m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);//ȫ����λ��Ŀ�����ֵ

															//m_GridPoint[i][j].Z = Well_Near_K(D);		//�ھ���Χ�ڵĿ�����ֵ20131212
															//m_GridPoint[i][j].Z = DisInv(D);			//�ھ���Χ�ڵľ��뷴�ȹ�ֵ20131212
				if (m_ZMin > m_GridPoint[i][j].Z)
					m_ZMin = m_GridPoint[i][j].Z;

				if (m_ZMax < m_GridPoint[i][j].Z)
					m_ZMax = m_GridPoint[i][j].Z;
			}
		}
		m_suV = Suv;		//20131212
		m_ni = ni;			//20130913 Ӧ�������Ʒ�
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
				//���ƴ˵������ֵ    û�жϲ�Ĺ���ֵ + ����ֵ
				//m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);//ȫ����λ��Ŀ�����ֵ
				//m_GridPoint[i][j].Z = Well_Near_K(D);		//�ھ���Χ�ڵĿ�����ֵ20131212
				m_GridPoint[i][j].Z = DisInv(D);			//�ھ���Χ�ڵľ��뷴�ȹ�ֵ20131212
				if (m_ZMin > m_GridPoint[i][j].Z)
					m_ZMin = m_GridPoint[i][j].Z;

				if (m_ZMax < m_GridPoint[i][j].Z)
					m_ZMax = m_GridPoint[i][j].Z;
			}

		}
	}

	//if ((int)m_oriData.size() <= 1)
	//{
	//	return ;		//����ֻ��һ����λ������20131107
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
	//		//���ƴ˵������ֵ    û�жϲ�Ĺ���ֵ + ����ֵ
	//		m_GridPoint[i][j].Z = InsertEst(Suv, D, ni);
	//		if(m_ZMin > m_GridPoint[i][j].Z)
	//			m_ZMin = m_GridPoint[i][j].Z;

	//		if(m_ZMax < m_GridPoint[i][j].Z)
	//			m_ZMax = m_GridPoint[i][j].Z;
	//	}

	//}

	/************************************************/
	//CString Str = "zֵ \n";
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
	//	str.Format("%f,    %f,    %f��   %f     ��%f\n",m_oriData[i].X,m_oriData[i].Y,m_oriData[i].Z,m_ZMax,m_ZMin);
	//	Str += str;
	//}
	//AfxMessageBox(Str);
	/************************************************/
	//m_suV = Suv;
	//m_ni = ni;				//20130913 Ӧ�������Ʒ�
}
void Jie_Ds::SetGridXY()
{
	//�����m_Border��Ӧ�ô洢���Ǿ��α߽�
	m_OriBoder = m_Border;	//ԭʼ�߽���Ϣ��m_OriBoder  20140806
	vector<Data> convexBag;
	//Charact();				//��ȡ�����Сֵ

							/*-------------------------------------------*/
							//���ݶ����߽粻�պϵ����20170718
	int OriCount = (int)m_OriBoder.size();
	if (OriCount >= 3)
	{
		Data Pt0 = m_OriBoder[0];
		Data Pt1 = m_OriBoder[OriCount - 1];
		if (GetDis(Pt0.X, Pt0.Y, Pt1.X, Pt1.Y) >= 0.0001)
		{
			//�ñ߽�պ�
			m_OriBoder.push_back(Pt0);
		}
	}

	/*-------------------------------------------*/

	/*-------------------------------------------*/
	//���ݾ�λ�������������m_B��ֵ 20150205  ���������ʲô����
	double k = sqrt(pow(m_XMax - m_XMin, 2) + pow(m_YMax - m_YMin, 2));
	if (m_B <= k)
	{
		m_B = k + 2;		//��Ҫ�Ƿ�ֹ�����Ǳ��ر��ķ�Χ
	}
	/*-------------------------------------------*/

	//GetRectBoder();		//ʹ�þ��α߽磨�Լ�����߽磩20131025

	//�ж��Ƿ��б߽磻����б߽磬�ʹӱ߽����͹������������ݵ����͹��
	//20131111m_Border��һֱ����ֵ�������ǵ����ľ�λ���ݻ����˹�����
	//�����ݵ����͹���ļ���
	if ((int)m_Border.size() <= 2)
	{
		m_Border.clear();
		//convexBag = Convex(m_oriData);
		convexBag = Withershins(m_oriData);	//�����µ�͹������20131108 ����͹��
		//ʹ��B��������͹���߽���Ż�
		//OptimizeBoder(convexBag,200);			
		//ʹ�ñ��������߽���͹���߽���Ż�
		OptimizeBorderBezier(convexBag,1000);
		m_Border = convexBag;
		m_OriBoder = m_Border;	//ԭʼ�߽���Ϣ��m_OriBoder  20140806
	}
	//�ɱ߽����͹���ļ���
	else
	{
		//cout << "hello2" << endl;
		convexBag = Withershins(m_Border);  //����20131107  ����͹��
		//ʹ��B��������͹���߽���Ż�
		//OptimizeBoder(convexBag);
		//ʹ�ñ��������߽���͹���߽���Ż�
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
	

	DividedFourParts((int)(convexBag.size() - 1), FourPort);//�Ż��߽��͹�����Ϊ������������

	//cout << "=========================================" << endl;
	//cout << "���͹��������һЩ��Ϣ" << endl;
	//cout << "FourPort[0]= " << FourPort[0] << " " << "FourPort[1]= " << FourPort[1] << " " << "FourPort[2]= " << FourPort[2] << " " << "FourPort[3]= " << FourPort[3] << endl;

	vector<Data> upLine, downLine, leftLine, rightLine;

	vector<Dis> upDis, downdis, leftDis, rightDis;

	//upLine��downline��leftLine��rightLine�洢��͹���Ĳ��ֱ߽��ϵĵ��λ��
	//upDis��downDis��leftDis��rightDis�е�dis�����洢�˶�Ӧ�㵽��һ����ľ���
	BordersPointDis(FourPort, convexBag, upLine, downLine, leftLine, rightLine,
		upDis, downdis, leftDis, rightDis);

	//����upLine��downline��leftLine��rightLine��Dis����ռ�ܵĳ��ȵİ�ֱ�
	//����upDis��downdis��leftDis��rightDis��angle��������x��������ļн�
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

	DL.clear();//������ٽ�����
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
//�õ�һ������������
double Jie_Ds::GetMagnitude(double fNumber)
{
	//������
	double magnitudeValue = 1.0;

	if (fNumber == 0.0)
		return(0.0);

	//�Ƿ�Ϊ����
	bool bNegative;
	bNegative = (fNumber<0) ? true : false;

	double positiveNumber = abs(fNumber);
	if (positiveNumber == 1)
	{//����1	
		magnitudeValue = 1.0;
	}
	else if (positiveNumber<1.0)
	{//С��1
		while (positiveNumber<1.0)
		{
			positiveNumber *= 10.0;
			magnitudeValue /= 10.0;
		}
	}
	else
	{//����1
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

	//�Ƿ�Ϊ����
	bool bNegative = false;
	bNegative = (StepMin<0) ? true : false;

	if (bNegative)
	{//�����ľ���ֵ�Ĺ���������֮���������෴
		bUporDown = !bUporDown;
	}

	//���Ȱ�����������
	StepMin = (float)fabs(StepMin);     //��С��0,��ȡ����ֵ

	if (!bUporDown)
	{//���¹���
	 //���㽥��ĸֵ
		dStep = (float)GetMagnitude(double(StepMin));//�õ��������������
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
	{//���Ϲ���
	 //���㽥��ĸֵ
		dStep = (float)GetMagnitude(double(StepMin));//�õ��������������
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

	//����������������ԭ�����
	if (bNegative)
	{
		RetVal *= -1;
	}
	return(RetVal);
}
void Jie_Ds::CalcSameArray()
{
	//��ֵ����С��ʾֵ
	m_Show_MinValue = FindStep((float)m_ZMin, false);//���¹���Сֵ

	//��ֵ�������ʾֵ
	m_Show_MaxValue = FindStep((float)(m_ZMax), true);//���Ϲ�����ֵ

																//��ֵ�߼��
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
//�õ����еĵ�ֵ��
//void Jie_Ds::EquivalentPoints(double Value, vector<THRVALUE>&Jie_VirtualIJK)
void Jie_Ds::EquivalentPoints(double Value, vector<TWOVALUE>&VirtualIJ)
{
	IsoLine temp_Line;
	VirtualIJ.clear();
	for (int i = 1; i <= m_XNum; i++)//�����ݱ��ϵĵ�ֵ��
	{
		for (int j = 1; j < m_YNum; j++)
		{
			//if (m_YFault[i][j])
			//{
			//	//���������н��� 20130913
			//	Y_FalutEquivalent(i, j, Value, VirtualIJ);//�õ�����ĵ�ֵ�� 20131015
			//	continue;
			//}
			//�жϴ˱����Ƿ��е�ֵ��
			if ((m_GridPoint[i][j].Z - Value) * (m_GridPoint[i][j + 1].Z - Value) < 0)
			{
				TWOVALUE P;
				P.X = i;
				//�����tֵһ�������ģ���������
				double t = ((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i][j + 1].Z - m_GridPoint[i][j].Z));
				P.Y = j + t;
				//if (t >= 1.0)
				//{
				//	int sk = 1;
				//	ASSERT(sk);
				//	CString str;
				//	str.Format("xt���� %f,%f,%f", Value, m_GridPoint[i][j].Z, m_GridPoint[i][j + 1].Z);
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
	for (int j = 1; j <= m_YNum; j++)//�����ϵĵ�ֵ��
	{
		for (int i = 1; i < m_XNum; i++)
		{
			//if (m_XFault[i][j])
			//{
			//	//���������н��� 20130913
			//	X_FalutEquivalent(i, j, Value, VirtualIJ);//�õ�����ĵ�ֵ�� 20131015
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
					//AfxMessageBox("yt����");
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

//����Ҫע�⣬����ͨ��vector<vector>����Ķ�ά���飬��һ��ֵ�����Ͻǣ��������½�
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
//ת��һ����
void Jie_Ds::GetReallyPoint(TWOVALUE A, TWOVALUE &B)
{
	double X1 = A.X;
	double Y1 = A.Y;
	int X0 = (int)X1;
	int Y0 = (int)Y1;
	double Mark = 0, C, D;
	//˵����һ�������
	if (X1 - X0 == 0 && Y1 - Y0 == 0)
	{
		C = m_GridPoint[X0][Y0].X;
		D = m_GridPoint[X0][Y0].Y;
	}
	//˵���ں���
	else if (X1 - X0 == 0 && Y1 - Y0 != 0)
	{
		Mark = abs(Y1 - Y0);
		C = (1 - Mark) * m_GridPoint[X0][Y0].X + Mark * (m_GridPoint[X0][Y0 + 1].X);//���������Բ�ֵ
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
	//if (m_NoFault[i][j]>-1)					//�������жϲ�
	//{
	//	TrackRightFault(Value, Line, IsoReal, X_Used, Y_Used);//׷��������
	//	return;
	//}
	//����׷��
	if ((A - Value) * (B - Value) < 0 && (D - Value) * (C - Value) < 0) //׷��������ֵ��
	{
		if (X_Used[i][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - A) / (A - B));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���±�
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j + 1] == 0 && ((D - Value) * (Z - Value) < 0 || ((B - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ϱ�
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}

	}
	else  //׷��һ����ֵ��
	{
		if ((A - Value) * (B - Value) < 0 && X_Used[i][j] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - A) / (A - B));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���±�
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if ((D - Value) * (C - Value) < 0 && X_Used[i][j + 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ϱ�
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && (B - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//����׷��
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
	//if (m_NoFault[i][j - 1]>-1)					//�������жϲ�
	//{
	//	TrackDownFault(Value, Line, IsoReal, X_Used, Y_Used);//׷��������
	//	return;
	//}
	//����׷��
	if ((A - Value) * (D - Value) < 0 && (B - Value) * (C - Value) < 0) //׷��������ֵ��
	{
		if (Y_Used[i][j - 1] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j - abs((Value - A) / (D - A));//�����൱��tmpV.Y = j -1+ abs((Value - D) / (A - D));��Ϊ-1+abs((Value - D) / (A - D))=-abs((Value - A) / (D - A))
			Line.Logic.push_back(tmpV);
			Y_Used[i][j - 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷�����
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
			//׷���ұ�
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else  //׷��һ����ֵ��
	{
		if ((A - Value) * (D - Value) < 0 && Y_Used[i][j - 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j - abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j - 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷�����
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
			//׷���ұ�
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j - 1] == 0 && (D - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j - 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j - 1] = 1;
			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);
			//����׷��
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
	//if (m_NoFault[i - 1][j]>-1)					//�������жϲ�
	//{
	//	TrackLeftFault(Value, Line, IsoReal, X_Used, Y_Used);//׷��������
	//	return;
	//}
	//����׷��
	if ((A - Value) * (B - Value) < 0 && (D - Value) * (C - Value) < 0) //׷��������
	{
		if (X_Used[i - 1][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - A) / (B - A));//�����൱��tmpV.X = i -1+ abs((Value - B) / (A - B));��Ϊ-1+abs((Value - B) / (A - B))=-abs((Value - A) / (B - A))
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���±�
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i - 1][j + 1] == 0 && ((B - Value) * (Z - Value) < 0 || ((D - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j + 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ϱ�
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else    //׷��һ����
	{
		if ((A - Value) * (B - Value) < 0 && X_Used[i - 1][j] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - A) / (B - A));
			tmpV.Y = j;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���±�
			TrackDown(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if ((D - Value) * (C - Value) < 0 && X_Used[i - 1][j + 1] == 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i - 1][j + 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ϱ�
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i - 1][j] == 0 && (B - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i - 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i - 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//����׷��
			//TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
			//�޸ģ�����׷��
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackUp(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{
	int Num = (int)Line.Logic.size();
	//��ǰ���һ����ֵ�����������
	int i = (int)Line.Logic[Num - 1].X;
	int j = (int)Line.Logic[Num - 1].Y;
	double A = m_GridPoint[i][j].Z;
	double B = m_GridPoint[i + 1][j].Z;
	double C = m_GridPoint[i + 1][j + 1].Z;
	double D = m_GridPoint[i][j + 1].Z;
	double Z = (A + B + C + D) / 4;
	//if (m_NoFault[i][j]>-1)					//�������жϲ�
	//{
	//	TrackUPFault(Value, Line, IsoReal, X_Used, Y_Used);//׷��������
	//	return;
	//}

	if ((A - Value) * (D - Value) < 0 && (B - Value) * (C - Value) < 0)//׷��������ֵ��
	{
		if (Y_Used[i][j] == 0 && ((A - Value) * (Z - Value) < 0 || ((C - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j + abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷�����
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && ((B - Value) * (Z - Value) < 0 || ((D - Value) * (Z - Value) < 0)))
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ұ�
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
	else	//ֻ��׷��һ����ֵ�� 
	{

		if (Y_Used[i][j] == 0 && (A - Value) * (D - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i;
			tmpV.Y = j + abs((Value - A) / (D - A));
			Line.Logic.push_back(tmpV);
			Y_Used[i][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷�����
			TrackLeft(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (Y_Used[i + 1][j] == 0 && (B - Value) * (C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + 1;
			tmpV.Y = j + abs((Value - B) / (C - B));
			Line.Logic.push_back(tmpV);
			Y_Used[i + 1][j] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//׷���ұ�
			TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		}
		else if (X_Used[i][j + 1] == 0 && (D - Value)*(C - Value) < 0)
		{
			TWOVALUE tmpV;
			tmpV.X = i + abs((Value - D) / (C - D));
			tmpV.Y = j + 1;
			Line.Logic.push_back(tmpV);
			X_Used[i][j + 1] = 1;

			//ʵ��ֵ
			TWOVALUE tmpR;
			GetReallyPoint(tmpV, tmpR);
			IsoReal.Logic.push_back(tmpR);

			//����׷��
			TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		}
	}
}
void Jie_Ds::TrackX(TWOVALUE A, double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used)
{

	int i = (int)A.X;
	int j = (int)A.Y;
	if (j > 1 && j < m_YNum)					//�Ǳ߽��,���ڵ�һ�������һ��֮��
	{
		//����ֵ
		IsoLine LineA;
		IsoLine LineB;
		LineA.Logic.push_back(A);
		LineB.Logic.push_back(A);


		IsoLine RLineA;
		IsoLine RLineB;


		//ʵ��ֵ
		TWOVALUE tmpR;
		//�õ���ʵ�ĵ�ֵ�������ֵ
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

		//TODO:����������ϵ����
		int NumA = (int)LineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //�������׵�
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

		//ʵ����������
		NumA = (int)RLineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA-1; k > 0; k--) //�������׵�
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
	else if (j == 1)//��һ��
	{

		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackUp(Value, Line, IsoReal, X_Used, Y_Used);
		X_Used[i][j] = 1;
	}
	else if (j == m_YNum)//���һ��
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
	if (i > 1 && i < m_XNum)					//�Ǳ߽�㣬���ڵ�һ�������һ��֮��
	{
		IsoLine LineA;
		IsoLine LineB;
		LineA.Logic.push_back(A);
		LineB.Logic.push_back(A);

		IsoLine RLineA;
		IsoLine RLineB;

		//ʵ��ֵ
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

		//������������
		int NumA = (int)LineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //�������׵�
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
		//ʵ����������
		NumA = (int)RLineA.Logic.size();
		if (NumA > 1)
		{
			for (int k = NumA - 1; k > 0; k--) //�������׵�
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
	else if (i == 1)//��һ��
	{
		Line.Logic.push_back(A);
		TWOVALUE tmpR;
		GetReallyPoint(A, tmpR);
		IsoReal.Logic.push_back(tmpR);
		TrackRight(Value, Line, IsoReal, X_Used, Y_Used);
		Y_Used[i][j] = 1;
	}
	else if (i == m_XNum)  //���һ��
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
//�ӳ�ʼֵ׷�ٵ�ֵ��
void Jie_Ds::TrackPoint(double Value, IsoLine &Line, IsoLine &IsoReal, vector<vector<double>>&X_Used, vector<vector<double>>&Y_Used, TWOVALUE First)
{
	//����֮����Ҫ������һ���жϣ�����Ϊ�����൱��һ��ȡ���Ĳ���
	//double t = abs((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i][j + 1].Z - m_GridPoint[i][j].Z));
	//double t = abs((Value - m_GridPoint[i][j].Z) / (m_GridPoint[i+1][j].Z - m_GridPoint[i][j].Z));
	//�����t�൱������������Բ�ֵǰ���ϵ����������֤t������һ��������
	int i = int(First.X);
	int j = (int)First.Y;
	//˵����ֵ���������
	if (i == First.X && j == First.Y)
	{
		//AfxMessageBox("��ֵ���������");			//���ڲ���
	}

	if (j == First.Y)								//����
	{
		if (X_Used[i][j] == 0)						//û�б�ʹ�ù�
		{
			//cout << "hello1" << endl;
			//TODO:׷�ٺ���
			//X_Used[int(First.X)][j] = 1;			//�׵㲻���жϣ������Ժ�׷�ٵ��׵�
			TrackX(First, Value, Line, IsoReal, X_Used, Y_Used);	//׷�ٺ��
		}
	}
	else											//����
	{
		if (Y_Used[i][j] == 0)						//û�б�ʹ�ù�
		{
			//cout << i << " " << j << endl;
			//cout << "hello2" << endl;
			//TODO:׷������
			//Y_Used[int(First.X)][j] = 1;			//�׵㲻���жϣ������Ժ�׷�ٵ��׵�
			TrackY(First, Value, Line, IsoReal, X_Used, Y_Used);	//׷���ݱ�
		}
	}
}
void Jie_Ds::TrackOneValue(double Value)
{
	vector<TWOVALUE> VirtualIJ;
	vector<THRVALUE> Jie_VirtualIJK;
	EquivalentPoints(Value, VirtualIJ);//�õ����еĵ�ֵ��
	vector<vector<double>> X_Used;
	vector<vector<double>> Y_Used;
    //����Ǹ�ֵ,����Ҫע������������ı�����ʽ
	SignBorder(X_Used, Y_Used);						
	//cout << "======================================================" << endl;
	//cout << "�����ǹ��ڱ������������Ϣ��" << endl;
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
	//��ֵ��׷��
	//cout << (int)VirtualIJ.size() << " ";
	for (int i = 0; i < (int)VirtualIJ.size(); i++)
	{
		IsoLine ISOLineOne;	//�洢��������
		IsoLine IsoReal;	//�洢ʵ������
		TWOVALUE First = VirtualIJ[i];
		TrackPoint(Value, ISOLineOne, IsoReal, X_Used, Y_Used, First);
		if (IsoReal.Logic.size() > 20)
		{
			ISOLineOne.Value = Value;
			IsoReal.Value = Value;
			m_IsoLine.push_back(ISOLineOne);
			m_IsoRealLine.push_back(IsoReal); //20131008���Զϲ��ڲ���������
		}
	}
	/*cout << "======================================================" << endl;
	cout << "������ִ�к���ڱ������������Ϣ��" << endl;
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

	//�ڴ˴���յ�ֵ��������
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
	Star.Y = 0;					//ƽ�ƣ����ж��Ƿ�һ��ֱ�߲����Ӱ�졤����������ΪС��~����߼��㾫��
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
//�߶��Ƿ��ཻ,index�������
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


	if (Xmax_1 < Xmin_2 || Xmin_1 > Xmax_2 || Ymin_1 > Ymax_2 || Ymax_1 < Ymin_2)   //���߶���С���β��ཻ���ó����߶β��ཻ
		return FALSE;
	else
	{
		V1 = (p1_e.X - p1_s.X) * (p2_s.Y - p1_s.Y) - (p2_s.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V2 = (p1_e.X - p1_s.X) * (p2_e.Y - p1_s.Y) - (p2_e.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V3 = (p2_e.X - p2_s.X) * (p1_s.Y - p2_s.Y) - (p1_s.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		V4 = (p2_e.X - p2_s.X) * (p1_e.Y - p2_s.Y) - (p1_e.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		//�����������ֱ����ʱ��   20131015
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
//��A�Ƿ���Line�ڲ�
bool Jie_Ds::IsInside(TWOVALUE A, vector<TWOVALUE> Line)
{
	//cout << "hhh" << endl;
	for (int i = 0; i < (int)Line.size() - 1; i++)
	{
		if (IsinLineK(Line[i], Line[i + 1], A))		//���ڱ߽���  20131023
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
	//�����ж��Ƿ����ڲ�
	TWOVALUE B = A;
	B.X = m_XMin - 200.0;     //AB����
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
	if (Sum % 2 == 1)//����Ϊ�����������ڲ�
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


	if (Xmax_1 < Xmin_2 || Xmin_1 > Xmax_2 || Ymin_1 > Ymax_2 || Ymax_1 < Ymin_2)   //���߶���С���β��ཻ���ó����߶β��ཻ
		return false;
	else   //���������Ĳ�����ʣ�������һ���߶ε������˵�����һ���߶ε�ͬһ��ʱ�����ཻ�������ཻ��
		   //( Q1 - P1 )��( P2 - P1) * ( P2 - P1)��(Q2 - P1) >= 0��
		   //����ļ��㹫ʽΪ:  P1 X P2 = x1y2 - x2y1;
	{
		V1 = (p1_e.X - p1_s.X) * (p2_s.Y - p1_s.Y) - (p2_s.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V2 = (p1_e.X - p1_s.X) * (p2_e.Y - p1_s.Y) - (p2_e.X - p1_s.X) * (p1_e.Y - p1_s.Y);
		V3 = (p2_e.X - p2_s.X) * (p1_s.Y - p2_s.Y) - (p1_s.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		V4 = (p2_e.X - p2_s.X) * (p1_e.Y - p2_s.Y) - (p1_e.X - p2_s.X) * (p2_e.Y - p2_s.Y);
		//�����������ֱ����ʱ��   20131015
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
	//cout << "��ǰ�����꣺" << OneIso.Logic[0].X << " " << OneIso.Logic[0].X << endl;
	if (OneIso.Logic.size() > 0)
	{
		bool flag = IsInside(OneIso.Logic[0], Jie_OriBoder);
		if (!flag)//�׵㲻��ԭʼ�߽��ڣ�ɾ����ǰ�ĵ㣬������߽��
		{
			for (int i = 1; i < (int)OneIso.Logic.size(); i++)
			{
				if (IsInside(OneIso.Logic[i], Jie_OriBoder))
				{
					TWOVALUE Star = OneIso.Logic[i - 1];
					TWOVALUE End = OneIso.Logic[i];
					TWOVALUE A;
					if (GetCrossPt(Star, End, Jie_OriBoder, A))//�õ��߽���������㹹�ɵ��߶εĽ���
					{
						TempLine.Logic.push_back(A);//�˽�����Ϊ��ֵ�ߵ����е�һ��
					}
					else
					{
						int sk = 0;
						//AfxMessageBox("111");
					}
					for (int j = i; j < (int)OneIso.Logic.size(); j++)
					{
						TempLine.Logic.push_back(OneIso.Logic[j]);//��������
					}
					break;
				}
			}
		}
		else
		{
			TempLine = OneIso;
		}
		if ((int)TempLine.Logic.size() <= 1)	//��ȫ����ԭʼ�߽�����
		{
			return;
		}
		IsoLine StarLine;
		StarLine.Value = TempLine.Value;
		StarLine.Logic.push_back(TempLine.Logic[0]);
		OneIso.Logic.clear();
		for (int i = 1; i < (int)TempLine.Logic.size(); i++)  //�ӵڶ����㿪ʼ��������
		{
			if (IsInside(TempLine.Logic[i], Jie_OriBoder))
			{
				StarLine.Logic.push_back(TempLine.Logic[i]);
			}
			else
			{
				//CString str;
				//str.Format("�׸��߽��%d,%f     %f",i,TempLine.Logic[0].X,TempLine.Logic[0].Y);
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
		if (OneIso.Logic.size() > 0)//������滹�����ݣ�˵����ֵ���м䲿��Ҳ���߽���з���
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
	Jie_IsoRealLine = NewIso;//���߽��и��ĵ�ֵ��
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