using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace SPA_AQ
{
    //public class OrbitalElements
    //{
    //    public double a;   // 半长轴
    //    public double e;   // 离心率
    //    public double i;   // 倾角
    //    public double raan;    // 升交点赤经
    //    public double w;   // 近地点幅角
    //    public double the; // 真近点角
    //    public double p;   // 半通径
    //    public double T;   // 轨道周期

    //    // 构造函数
    //    public OrbitalElements(double a, double e, double i, double raan, double w, double the)
    //    {
    //        this.a = a;
    //        this.e = e;
    //        this.i = i;
    //        this.raan = raan;
    //        this.w = w;
    //        this.the = the;
    //        this.p = this.a * (1 - this.e * this.e);
    //        this.T = GetT(398600.5e9);
    //    }

    //    // 计算轨道周期
    //    public double GetT(double mu)
    //    {
    //        return 2 * Math.PI * Math.Sqrt(a * a * a / mu);
    //    }

    //    // 获得轨道参数向量
    //    public double[] GetElement0()
    //    {
    //        return new double[] { a, e, i, raan, w, the };
    //    }
    //}

    // 工具类
    public class TwoBody
    {
        // 惯性系 to 轨道系
        // 输入参数： orbital elemet OrbitalElements
        // 输出参数：位置矢量 R0、速度矢量 V0
        //public static void Orb2Inert(OrbitalElements elements, out Vector<double> R0, out Vector<double> V0)
        //{
        //    //轨道根数到地心惯性系
        //    double miu = 398600.5 * Math.Pow(10, 9);

        //    double r0 = elements.p / (1 + elements.e * Math.Cos(elements.the));
        //    Vector<double> FR0 = Vector<double>.Build.Dense(new[] { r0 * Math.Cos(elements.the), r0 * Math.Sin(elements.the), 0 });

        //    double vP = Math.Sqrt(miu / (elements.p * (1 - Math.Pow(elements.e, 2))));
        //    Vector<double> FV0 = Vector<double>.Build.Dense(new[] { -vP * Math.Sin(elements.the), vP * (Math.Cos(elements.the) + elements.e), 0 });

        //    Matrix<double> Cab = Matrix<double>.Build.DenseOfArray(new[,]
        //    {
        //    { Math.Cos(elements.raan) * Math.Cos(elements.w) - Math.Sin(elements.raan) * Math.Sin(elements.w) * Math.Cos(elements.i), -Math.Cos(elements.raan) * Math.Sin(elements.w) - Math.Sin(elements.raan) * Math.Cos(elements.w) * Math.Cos(elements.i), Math.Sin(elements.raan) * Math.Sin(elements.i) },
        //    { Math.Sin(elements.raan) * Math.Cos(elements.w) + Math.Cos(elements.raan) * Math.Sin(elements.w) * Math.Cos(elements.i), -Math.Sin(elements.raan) * Math.Sin(elements.w) + Math.Cos(elements.raan) * Math.Cos(elements.w) * Math.Cos(elements.i), -Math.Cos(elements.raan) * Math.Sin(elements.i) },
        //    { Math.Sin(elements.w) * Math.Sin(elements.i), Math.Cos(elements.w) * Math.Sin(elements.i), Math.Cos(elements.i) }
        //});

        //    R0 = Cab.Multiply(FR0);
        //    V0 = Cab.Multiply(FV0);
        //}

        public static void ElementsToRV(double a, double e, double i, double RAAN, double w, double theta, out Vector<double> R, out Vector<double> V)
        {   //功能其实和Orb2Inert一样，但是历史遗留问题导致，主程序中一部分输入为结构体，一部分为分量

            double mu = 398600.5 * Math.Pow(10, 9);
            double p = a * (1 - Math.Pow(e, 2));
            double r0 = p / (1 + e * Math.Cos(theta));
            Vector<double> FR0 = Vector<double>.Build.Dense(new[] { r0 * Math.Cos(theta), r0 * Math.Sin(theta), 0 });
            double vP = Math.Sqrt(mu / p);
            Vector<double> FV0 = Vector<double>.Build.Dense(new[] { -vP * Math.Sin(theta), vP * (Math.Cos(theta) + e), 0 });

            Matrix<double> Cab = Matrix<double>.Build.DenseOfArray(new[,]
            {
            { Math.Cos(RAAN) * Math.Cos(w) - Math.Sin(RAAN) * Math.Sin(w) * Math.Cos(i), -Math.Cos(RAAN) * Math.Sin(w) - Math.Sin(RAAN) * Math.Cos(w) * Math.Cos(i), Math.Sin(RAAN) * Math.Sin(i) },
            { Math.Sin(RAAN) * Math.Cos(w) + Math.Cos(RAAN) * Math.Sin(w) * Math.Cos(i), -Math.Sin(RAAN) * Math.Sin(w) + Math.Cos(RAAN) * Math.Cos(w) * Math.Cos(i), -Math.Cos(RAAN) * Math.Sin(i) },
            { Math.Sin(w) * Math.Sin(i), Math.Cos(w) * Math.Sin(i), Math.Cos(i) }
        });
            //近焦点系in to cartesian
            R = Cab.Multiply(FR0);
            V = Cab.Multiply(FV0);
        }

        public static void InertToOrbR(Vector<double> RI, double aT, double eT, double iT, double RAANT, double wT, double thetaT,
                                        out Vector<double> Ro)
        {
            //位置转到（第二）轨道坐标系中描述,其中带T的指建立轨道系的卫星的六根数
            double mu = 398600.5 * Math.Pow(10, 9);
            double re = 6371.11 * Math.Pow(10, 3);

            // 求惯性系到轨道系坐标转换矩阵
            double uT = wT + thetaT;
            Matrix<double> COA = Matrix<double>.Build.DenseOfArray(new[,]
            {
            { -Math.Sin(uT) * Math.Cos(RAANT) - Math.Cos(uT) * Math.Cos(iT) * Math.Sin(RAANT), -Math.Sin(uT) * Math.Sin(RAANT) + Math.Cos(uT) * Math.Cos(iT) * Math.Cos(RAANT), Math.Cos(uT) * Math.Sin(iT) },
            { -Math.Sin(iT) * Math.Sin(RAANT), Math.Sin(iT) * Math.Cos(RAANT), -Math.Cos(iT) },
            { -Math.Cos(uT) * Math.Cos(RAANT) + Math.Sin(uT) * Math.Cos(iT) * Math.Sin(RAANT), -Math.Cos(uT) * Math.Sin(RAANT) - Math.Sin(uT) * Math.Cos(iT) * Math.Cos(RAANT), -Math.Sin(uT) * Math.Sin(iT) }

            });
            // 求惯性系到轨道系坐标转换矩阵结束
            //坐标转换
            Ro = COA * RI;

        }

        public static void PositionSunOrb(double t, Vector<double> ET0, Vector<double> RS_I, out Vector<double> R_o)
        {   //求解t时刻太阳在目标轨道系下的位置 （m）
            //% ET0为t0时刻目标的 orbital elemet
            //% RS_I为t0时刻太阳在地心J2000坐标系下的位置(m)
            //% t0时刻为2025年4.18 9:00
            double mu = 398600.5 * Math.Pow(10, 9);

            // t 时刻目标的 orbital elemet
            double aT0 = ET0[0];
            double nT = Math.Sqrt(mu / Math.Pow(aT0, 3));
            double thetaT = ET0[5] + nT * t;

            // t 时刻目标的位置速度
            TwoBody.ElementsToRV(ET0[0], ET0[1], ET0[2], ET0[3], ET0[4], thetaT, out Vector<double> RT, out Vector<double> VT); // t 时刻目标的位置速度

            Vector<double> R_T2S = RS_I - RT; // J2000 下目标指向太阳位置向量
            TwoBody.InertToOrbR(R_T2S, ET0[0], ET0[1], ET0[2], ET0[3], ET0[4], thetaT, out R_o); // 目标轨道系下目标-太阳位置矢量描述
        }

        public static void OrbToInertRV(Vector<double> Ro, Vector<double> Vo, double aT, double eT, double iT, double RAANT, double wT, double thetaT,
                                        out Vector<double> RI, out Vector<double> VI)
        {
            // 将在第二轨道坐标系下描述的位置、速度转到惯性系中描述，其中带T的指建立轨道系的卫星的六根数

            double mu = 398600.5 * Math.Pow(10, 9);
            // 求定义轨道系卫星的位置速度矢量
            double nT = Math.Sqrt(mu / Math.Pow(aT, 3));
            double pT = aT * (1 - Math.Pow(eT, 2));
            double rT = pT / (1 + eT * Math.Cos(thetaT));
            Vector<double> FRT = Vector<double>.Build.DenseOfArray(new[] { rT * Math.Cos(thetaT), rT * Math.Sin(thetaT), 0 });
            Vector<double> FVT = Vector<double>.Build.DenseOfArray(new[] { -Math.Sqrt(mu / pT) * Math.Sin(thetaT), Math.Sqrt(mu / pT) * (Math.Cos(thetaT) + eT), 0 });

            //Cab为近焦点坐标系in to cartesian
            Matrix<double> CabT = Matrix<double>.Build.DenseOfArray(new[,]
            {
            { Math.Cos(RAANT) * Math.Cos(wT) - Math.Sin(RAANT) * Math.Sin(wT) * Math.Cos(iT), -Math.Cos(RAANT) * Math.Sin(wT) - Math.Sin(RAANT) * Math.Cos(wT) * Math.Cos(iT), Math.Sin(RAANT) * Math.Sin(iT) },
            { Math.Sin(RAANT) * Math.Cos(wT) + Math.Cos(RAANT) * Math.Sin(wT) * Math.Cos(iT), -Math.Sin(RAANT) * Math.Sin(wT) + Math.Cos(RAANT) * Math.Cos(wT) * Math.Cos(iT), -Math.Cos(RAANT) * Math.Sin(iT) },
            { Math.Sin(wT) * Math.Sin(iT), Math.Cos(wT) * Math.Sin(iT), Math.Cos(iT) }
        });

            Vector<double> RT = CabT * FRT;
            Vector<double> VT = CabT * FVT;
            //求定义轨道系卫星的位置速度矢量结束


            // 惯性系下的位置速度
            double uT = wT + thetaT;
            Matrix<double> Coa = Matrix<double>.Build.DenseOfArray(new[,]
            {
            { -Math.Sin(uT) * Math.Cos(RAANT) - Math.Cos(uT) * Math.Cos(iT) * Math.Sin(RAANT), -Math.Sin(uT) * Math.Sin(RAANT) + Math.Cos(uT) * Math.Cos(iT) * Math.Cos(RAANT), Math.Cos(uT) * Math.Sin(iT) },
            { -Math.Sin(iT) * Math.Sin(RAANT), Math.Sin(iT) * Math.Cos(RAANT), -Math.Cos(iT) },
            { -Math.Cos(uT) * Math.Cos(RAANT) + Math.Sin(uT) * Math.Cos(iT) * Math.Sin(RAANT), -Math.Cos(uT) * Math.Sin(RAANT) - Math.Sin(uT) * Math.Cos(iT) * Math.Cos(RAANT), -Math.Sin(uT) * Math.Sin(iT) }
        });

            RI = RT + Coa.Transpose() * Ro;
            VI = VT + Coa.Transpose() * (Vo + Vector<double>.Build.DenseOfArray(new[] { -nT * Ro[2], 0, nT * Ro[0] }));

        }

        public static void PositionInInertial(double t, Vector<double> ET0, Vector<double> Ro, Vector<double> Vo,
                                       out Vector<double> Rf, out Vector<double> Vf)
        {
            //由时刻t(s)确定惯性系下的终端位置速度(m,m/s)
            //ET0为t0时刻目标的 orbital elemet
            double mu = 398600.5 * Math.Pow(10, 9);
            double aT0 = ET0[0];
            double nT = Math.Sqrt(mu / Math.Pow(aT0, 3));
            double thetaT = ET0[5] + nT * t;

            OrbToInertRV(Ro, Vo, ET0[0], ET0[1], ET0[2], ET0[3], ET0[4], thetaT, out Rf, out Vf);
        }

        


        public static void RelativeLambert(double n, double t1, double t2, Vector<double> R1, Vector<double> R2, out Vector<double> V1, out Vector<double> V2)
        {   //求解相对Lambert问题
            double t0 = 0.0;
            double t = t2 - t1 + t0;
            double tau = n * (t - t0);

            Matrix<double> Phirr = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 1.0, 0.0, 6.0 * (tau - Math.Sin(tau)) },
                { 0.0, Math.Cos(tau), 0 },
                { 0.0, 0.0, 4.0 - 3.0 * Math.Cos(tau) }
            });

            Matrix<double> Phirv = Matrix<double>.Build.DenseOfArray(new double[,] {
                { (4.0 * Math.Sin(tau) - 3.0 * tau) / n, 0.0, 2.0 * (1.0 - Math.Cos(tau)) / n },
                { 0.0, Math.Sin(tau) / n, 0.0 },
                { 2.0 * (Math.Cos(tau) - 1.0) / n, 0.0, Math.Sin(tau) / n }
            });

            Matrix<double> Phivr = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 0.0, 0.0, 6.0 * n * (1.0 - Math.Cos(tau)) },
                { 0.0, -n * Math.Sin(tau), 0.0 },
                { 0.0, 0.0, 3.0 * n * Math.Sin(tau) }
            });

            Matrix<double> Phivv = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 4.0 * Math.Cos(tau) - 3.0, 0.0, 2.0 * Math.Sin(tau) },
                { 0.0, Math.Cos(tau), 0.0 },
                { -2.0 * Math.Sin(tau), 0.0, Math.Cos(tau) }
            });

            Matrix<double> Grr = -Phirv.Inverse() * Phirr;
            Matrix<double> Grv = Phirv.Inverse();
            Matrix<double> Gvr = Phivr - Phivv * Grv * Phirr;
            Matrix<double> Gvv = Phivv * Grv;

            V1 = Grr * R1 + Grv * R2;
            V2 = Gvr * R1 + Gvv * R2;
        }

        // 测试 X_cw_clfm_propagate
        //static void Main(string[] args)
        //{
        //    // Test inputs
        //    double t = 72600;
        //    double[] tp = Enumerable.Range(0, 13).Select(i => 72000.0 + i * 7200.0).ToArray();
        //    double t0 = 0.0;
        //    double t_app = 72000.0;
        //    Vector<double> X0 = Vector<double>.Build.DenseOfArray(new double[] { 6.853241613188095e+04, -1.316882960073963e+04, -5.467162524874004e+03, 0.0, 0.0, 0.0 });
        //    Matrix<double> DeltaV = Matrix<double>.Build.DenseOfArray(new double[,] {
        //        { 0.721568742877672, 1.93287262711887, 1.00278383950403, -0.197355155665785, -1.34434732186935, -2.12929974063193, -2.34077370843739, -1.92174789872202, -0.985038402693204, 0.217039943366320, 1.36064682199899, 2.13770286159094 },
        //        { -0.258269730686071, -0.516011619474367, -0.516011539981945, -0.516011511825270, -0.516011563329312, -0.516011642686583, -0.516011670123616, -0.516011618124074, -0.516011538970475, -0.516011512146432, -0.516011564575869, -0.516011643629207 },
        //        { -4.18128609950665, 2.67122982463454, 4.24985086086432, 4.68304249784307, 3.85388366480346, 1.98614893636015, -0.416408536972736, -2.70642507557813, -4.26752850300233, -4.67991109037711, -3.83272645938034, -1.95371642703859 }
        //    });
        //    double nT = 7.291641582282420e-05;
        //    int N = 12;

        //    // Call X_cw_clfm_propagate method
        //    Vector<double> X = X_cw_clfm_propagate(t, tp, t0, t_app, X0, DeltaV, nT, N);

        //    // Print results
        //    Console.WriteLine(X.ToString());
        //    Console.ReadKey();

        //}


        // 测试程序 InertToOrbR
        //static void Main(string[] args)
        //{
        //    Vector<double> RI = Vector<double>.Build.Dense(new[] { 132323868556.545, 65228226287.6864, 28258311697.5917 });
        //    double aT = 42166000;
        //    double eT = 0;
        //    double iT = 0;
        //    double RAANT = 0;
        //    double wT = 0;
        //    double thetaT = 5.249981939243343;
        //    TwoBody.InertToOrbR(RI, aT, eT, iT, RAANT, wT, thetaT,out Vector<double> Ro);

        //    Console.WriteLine($"Ro_x = {Ro[0]}");
        //    Console.WriteLine($"Ro_y = {Ro[1]}");
        //    Console.WriteLine($"Ro_z = {Ro[2]}");

        //    Console.ReadKey();
        //}

        // 测试PositionSunOrb
        //static void Main(string[] args)
        //{
        //    double t = 72000;
        //    Vector<double> ET0 = Vector<double>.Build.Dense(new[] { 42166000, 0, 0.0, 0.0, 0.0, 0.0 });
        //    Vector<double> RS_I = Vector<double>.Build.Dense(new[] { 1.323454604986535e+11, 6.519200806586389e+10, 2.825831169759174e+10 });
        //    TwoBody.PositionSunOrb(t, ET0, RS_I,out Vector<double> RS_o);
        //    Console.WriteLine($"Ro_x = {RS_o[0]}");
        //    Console.WriteLine($"Ro_y = {RS_o[1]}");
        //    Console.WriteLine($"Ro_z = {RS_o[2]}");

        //    Console.ReadKey();
        //}

        // 测试PositionInInertial
        //    static void Main(string[] args)
        //{
        //    double t = 72000;
        //    Vector<double> ET0 = Vector<double>.Build.Dense(new[] { 42166000, 0, 0.0, 0.0, 0.0, 0.0 });
        //    Vector<double> Ro = Vector<double>.Build.Dense(new[] { 68532.4161318810, -13168.8296007396, -5467.16252487400 });
        //    Vector<double> Vo = Vector<double>.Build.Dense(new[] { 0,0.0,0.0 });

        //    TwoBody.PositionInInertial(t, ET0, Ro, Vo, out Vector<double> Rf, out Vector<double> Vf);
        //    Console.WriteLine($"Rf_x = {Rf[0]}");
        //    Console.WriteLine($"Rf_y = {Rf[1]}");
        //    Console.WriteLine($"Rf_z = {Rf[2]}");
        //    Console.WriteLine($"Vf_x = {Vf[0]}");
        //    Console.WriteLine($"Vf_y = {Vf[1]}");
        //    Console.WriteLine($"Vf_z = {Vf[2]}");
        //    Console.ReadKey();
        //}

        //测试 SolveLambert
        //public static void Main(string[] args)
        //{
        //    Vector<double> R1 = Vector<double>.Build.Dense(new[] {4.177391587524976e+07, 2.133180430980697e+06,-2.382939839144658e+05});
        //    Vector<double> R2 = Vector<double>.Build.Dense(new[] { 42087103.1484871, 2578236.74770519, 0 });
        //    Vector<double> V1b = Vector<double>.Build.Dense(new[] { -1269.60096116431, 4026.67545536416, -768.312414226506 });
        //    Vector<double> V2a = Vector<double>.Build.Dense(new[] { -187.995780576534, 3068.84071395083, 0 });
        //    double t1 = 7.602824215528237e+02;
        //    double t2 = 8.553475192120985e+02;

        //    SolveLambertMultiRevolution(R1, R2, t1, t2, V1b, V2a, out double flag, out Vector<double> V1, out Vector<double> V2);

        //    Console.WriteLine($"rf_x = { V1[0]}");
        //    Console.WriteLine($"rf_y = { V1[1]}");
        //    Console.WriteLine($"rf_z = { V1[2]}");
        //    Console.WriteLine($"vf_x = {V2[0]}");
        //    Console.WriteLine($"vf_y = {V2[1]}");
        //    Console.WriteLine($"vf_z = {V2[2]}");
        //    Console.ReadKey();
        //}

        // 绕飞测试
        //static void Main(string[] args)
        //{
        //    Vector<double> RS_I = Vector<double>.Build.Dense(new[] { 1.323454604986535e+11, 6.519200806586389e+10, 2.825831169759174e+10 }); // 太阳在地心惯性系下位置
        //    double miu = 398600.5 * Math.Pow(10, 9);   // 地心引力常数
        //    double re = 6371.11 * Math.Pow(10, 3);     // 地球半径
        //    double delta = 5 * Math.PI / 180;          // 半顶角
        //    double T_fly = 24 * 3600.0;                 // 保持时间
        //    int N = 12;                                // 保持次数
        //    double t0 = 0;

        //    // space shuttle orbital elemet
        //    double aP = 42.166 * Math.Pow(10, 6), eP = 0, iP = 0, RAANP = 0, wP = 0, thetaP = -0.6794 * Math.PI / 180;
        //    double pP = aP * (1 - eP * eP);
        //    double TP = 2 * Math.PI * Math.Sqrt(aP * aP * aP / miu);
        //    double[] et0 = new double[] { aP, eP, iP, RAANP, wP, thetaP };
        //    //space shuttle orbital elemetin to cartesianR0 V0
        //    ElementsToRV(aP, eP, iP, RAANP, wP, thetaP, out Vector<double> RP0, out Vector<double> VP0);

        //    //ISS  orbital elemet
        //    double aT = 42.166 * Math.Pow(10, 6), eT = 0, iT = 0, RAANT = 0, wT = 0, thetaT = 0;
        //    double pT = aT * (1 - eT * eT);
        //    double TT = 2 * Math.PI * Math.Sqrt(aT * aT * aT / miu);

        //    Vector<double> ET0 = Vector<double>.Build.Dense(new[] { aT, eT, iT, RAANT, wT, thetaT });

        //    // LONG DIST APPR
        //    double t_app = 20 * 3600;
        //    double r_fly = 70e3;
        //    PositionSunOrb(t_app, ET0, RS_I, out Vector<double> RS_o);

        //    var Rapp_o = RS_o.Normalize(2) * r_fly;
        //    Vector<double> Vapp_o = Vector<double>.Build.Dense(new[] { 0.0, 0, 0 });
        //    // IN TO  icrf
        //    PositionInInertial(t_app, ET0, Rapp_o, Vapp_o,
        //                               out Vector<double> R_app, out Vector<double> V_app);
        //    // 远距离Lambert抵近
        //    SolveLambert(RP0, R_app, t0, t_app, VP0,
        //                                 out double flg, out Vector<double> VPa, out Vector<double> Vfb);
        //    //接近段的消耗
        //    double dv_app = (VPa - VP0).L2Norm() + (V_app - Vfb).L2Norm();
        //    // 受迫绕飞
        //    //ISS轨道角速率(rad/s)
        //    double nT = Math.Sqrt(miu / (aP * aP * aP));
        //    //中间时刻与中间点
        //    double[] tp = new double[N + 1];
        //    Matrix<double> Rp = Matrix<double>.Build.Dense(3, N + 1);

        //    tp[0] = t_app;
        //    Rp.SetColumn(0, Rapp_o);

        //    for (int i = 1; i <= N; i++)
        //    {
        //        tp[i] = t_app + i * T_fly / N;
        //        PositionSunOrb(tp[i], ET0, RS_I, out Vector<double> RS_o_i);
        //        Rp.SetColumn(i, RS_o_i.Normalize(2) * r_fly);
        //    }
        //    //Lambert转移
        //    //每一列代表每段Lambert转移的初始速度
        //    Matrix<double> Vpa = Matrix<double>.Build.Dense(3, N);
        //    //每一列代表每段Lambert转移的末端速度
        //    Matrix<double> Vpb = Matrix<double>.Build.Dense(3, N);
        //    for (int i = 0; i < N; i++)
        //    {
        //        RelativeLambert(nT, tp[i], tp[i + 1], Rp.Column(i), Rp.Column(i + 1),
        //                    out Vector<double> V1, out Vector<double> V2);
        //        Vpa.SetColumn(i, V1);
        //        Vpb.SetColumn(i, V2);
        //    }
        //    //整个过程的消耗及轨迹
        //    // 脉冲
        //    Matrix<double> DeltaV = Matrix<double>.Build.Dense(3, N);
        //    double[] deltav = new double[N];
        //    DeltaV.SetColumn(0, Vpa.Column(0) - Vapp_o);
        //    deltav[0] = DeltaV.Column(0).L2Norm();
        //    for (int i = 1; i < N; i++)
        //    {
        //        DeltaV.SetColumn(i, Vpa.Column(i) - Vpb.Column(i - 1));
        //        deltav[i] = DeltaV.Column(i).L2Norm();
        //    }
        //    double sdv_fa = deltav.Sum();
        //    double sumdv = dv_app + sdv_fa;
        //    //Console.WriteLine($"sumdv = {sumdv}");
        //    //Console.ReadKey();

        //    //变化趋势
        //    // 生成目标的位置序列
        //    double[] t_target_array = new double[2 * N + 1];
        //    for (int i = 0; i <= 2 * N; i++)
        //    {
        //        double t = t_app + i * T_fly / (2 * N);
        //        t_target_array[i] = t;
        //    }

        //    Matrix<double> L_TS_v = Matrix<double>.Build.Dense(3, t_target_array.Length);

        //    for (int i = 0; i < t_target_array.Length; i++)
        //    {
        //        PositionSunOrb(t_target_array[i], ET0, RS_I, out Vector<double> RS_o_t);
        //        L_TS_v.SetColumn(i, RS_o_t.Normalize(2));
        //    }

        //    //变化图象绘制，这里开始准备space shuttle的绘图变量，包括space shuttle的位置序列
        //    double[] t_array = new double[(int)Math.Ceiling(T_fly / 120) + 1];
        //    for (int i = 0; i < t_array.Length; i++)
        //    {
        //        t_array[i] = t_app + i * 120;
        //    }
        //    Matrix<double> X_cw_array = Matrix<double>.Build.Dense(6, t_array.Length);
        //    Matrix<double> R = Matrix<double>.Build.Dense(3, t_array.Length);
        //    double[] cos_the = new double[t_array.Length];
        //    double[] Fr = new double[t_array.Length];

        //    Vector<double> X0 = Vector<double>.Build.DenseOfArray(new double[] { Rapp_o[0], Rapp_o[1], Rapp_o[2], Vapp_o[0], Vapp_o[1], Vapp_o[2] });

        //    for (int i = 0; i < t_array.Length; i++)
        //    {
        //        double t = t_array[i];
        //        Vector<double> X_cw = X_cw_clfm_propagate(t, tp, t0, t_app, X0, DeltaV, nT, N);

        //        X_cw_array.SetColumn(i, X_cw);
        //        Vector<double> R_vector = (Matrix<double>.Build.DenseOfArray(new double[,]
        //        {
        //            {1, 0, 0, 0, 0, 0},
        //            {0, 1, 0, 0, 0, 0},
        //            {0, 0, 1, 0, 0, 0}
        //        }) * X_cw).SubVector(0, 3);
        //        R.SetColumn(i, R_vector);

        //        Console.WriteLine($"sumdv = {R.Column(i)}");
        //        Console.ReadKey();
        //    }



        //}
    }
}
