using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;


namespace SPA_AQ
{
    class frontlight
    {
        public static Vector<double> X_cw_clfm_propagate(double t, double[] tp, double t0, double t_app, Vector<double> X0, Matrix<double> DeltaV, double nT, int N)
        {
            double tau = nT * (t - t_app + t0 - t0);

            Matrix<double> F_rr = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 1.0, 0.0, 6.0 * (tau - Math.Sin(tau)) },
                { 0.0, Math.Cos(tau), 0.0 },
                { 0.0, 0.0, 4.0 - 3.0 * Math.Cos(tau) }
            });

            Matrix<double> F_rv = Matrix<double>.Build.DenseOfArray(new double[,] {
                { (4.0 * Math.Sin(tau) - 3.0 * tau) / nT, 0.0, 2.0 * (1.0 - Math.Cos(tau)) / nT },
                { 0.0, Math.Sin(tau) / nT, 0.0 },
                { 2.0 * (Math.Cos(tau) - 1.0) / nT, 0.0, Math.Sin(tau) / nT }
            });

            Matrix<double> F_vr = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 0.0, 0.0, 6.0 * nT * (1.0 - Math.Cos(tau)) },
                { 0.0, -nT * Math.Sin(tau), 0.0 },
                { 0.0, 0.0, 3.0 * nT * Math.Sin(tau) }
            });

            Matrix<double> F_vv = Matrix<double>.Build.DenseOfArray(new double[,] {
                { 4.0 * Math.Cos(tau) - 3.0, 0.0, 2.0 * Math.Sin(tau) },
                { 0.0, Math.Cos(tau), 0.0 },
                { -2.0 * Math.Sin(tau), 0.0, Math.Cos(tau) }
            });

            Matrix<double> F = Matrix<double>.Build.Dense(6, 6);
            F.SetSubMatrix(0, 0, F_rr);
            F.SetSubMatrix(0, 3, F_rv);
            F.SetSubMatrix(3, 0, F_vr);
            F.SetSubMatrix(3, 3, F_vv);

            Vector<double> X = F * X0;

            for (int i = 0; i < N; i++)
            {
                double tauv = nT * (t - tp[i] + t0 - t0);

                Matrix<double> F_rv_i = Matrix<double>.Build.DenseOfArray(new double[,] {
            { (4.0 * Math.Sin(tauv) - 3.0 * tauv) / nT, 0.0, 2.0 * (1.0 - Math.Cos(tauv)) / nT },
            { 0.0, Math.Sin(tauv) / nT, 0.0 },
            { 2.0 * (Math.Cos(tauv) - 1.0) / nT, 0.0, Math.Sin(tauv) / nT }
        });

                Matrix<double> F_vv_i = Matrix<double>.Build.DenseOfArray(new double[,] {
            { 4.0 * Math.Cos(tauv) - 3.0, 0.0, 2.0 * Math.Sin(tauv) },
            { 0.0, Math.Cos(tauv), 0.0 },
            { -2.0 * Math.Sin(tauv), 0.0, Math.Cos(tauv) }
        });

                Matrix<double> F_V = F_rv_i.Transpose().Append(F_vv_i.Transpose()).Transpose();

                X += F_V * DeltaV.Column(i) * (t >= tp[i] ? 1.0 : 0.0);
            }
            return X;
        }

        //要使用顺光时，将这里注释取消
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
        //    TwoBody.ElementsToRV(aP, eP, iP, RAANP, wP, thetaP, out Vector<double> RP0, out Vector<double> VP0);

        //    //ISS  orbital elemet
        //    double aT = 42.166 * Math.Pow(10, 6), eT = 0, iT = 0, RAANT = 0, wT = 0, thetaT = 0;
        //    double pT = aT * (1 - eT * eT);
        //    double TT = 2 * Math.PI * Math.Sqrt(aT * aT * aT / miu);

        //    Vector<double> ET0 = Vector<double>.Build.Dense(new[] { aT, eT, iT, RAANT, wT, thetaT });

        //    // LONG DIST APPR
        //    double t_app = 20 * 3600;
        //    double r_fly = 70e3;
        //    TwoBody.PositionSunOrb(t_app, ET0, RS_I, out Vector<double> RS_o);

        //    var Rapp_o = RS_o.Normalize(2) * r_fly;
        //    Vector<double> Vapp_o = Vector<double>.Build.Dense(new[] { 0.0, 0, 0 });
        //    // IN TO  icrf
        //    TwoBody.PositionInInertial(t_app, ET0, Rapp_o, Vapp_o,
        //                               out Vector<double> R_app, out Vector<double> V_app);
        //    // 远距离Lambert抵近
        //    TwoBody.SolveLambert(RP0, R_app, t0, t_app, VP0,
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
        //        TwoBody.PositionSunOrb(tp[i], ET0, RS_I, out Vector<double> RS_o_i);
        //        Rp.SetColumn(i, RS_o_i.Normalize(2) * r_fly);
        //    }
        //    //Lambert转移
        //    //每一列代表每段Lambert转移的初始速度
        //    Matrix<double> Vpa = Matrix<double>.Build.Dense(3, N);
        //    //每一列代表每段Lambert转移的末端速度
        //    Matrix<double> Vpb = Matrix<double>.Build.Dense(3, N);
        //    for (int i = 0; i < N; i++)
        //    {
        //        TwoBody.RelativeLambert(nT, tp[i], tp[i + 1], Rp.Column(i), Rp.Column(i + 1),
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
        //        TwoBody.PositionSunOrb(t_target_array[i], ET0, RS_I, out Vector<double> RS_o_t);
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
