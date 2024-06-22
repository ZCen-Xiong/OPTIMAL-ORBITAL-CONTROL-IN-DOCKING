using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace SPA_AQ
{
    class Multi
    {
        public static void CombineSolvePropagation(double t1, Vector<double> R0, Vector<double> V0, double t0, out Vector<double> R1, out Vector<double> V1)
        {
            // 普适变量+状态转移函数轨道预报
            // 用二分法+牛顿下山法解普适变量x
            double mu = 398600.5 * Math.Pow(10, 9);

            if (Math.Abs(t1 - t0) <= Math.Pow(10, -4))
            {
                R1 = R0;
                V1 = V0;
            }
            else
            {
                double r0 = R0.Norm(2);
                // vector类型没有叉乘，只有vector3有，但是后面用到不少vector的类，考虑到叉乘只在这个子函数中用了，所以直接硬写
                var R0Vec3 = new Vector3((float)R0[0], (float)R0[1], (float)R0[2]);
                var V0Vec3 = new Vector3((float)V0[0], (float)V0[1], (float)V0[2]);   //vector3只能float
                var R0XV0_f = Vector3.Cross(R0Vec3, V0Vec3);
                var R0XV0Array = new double[3] { R0XV0_f.X, R0XV0_f.Y, R0XV0_f.Z };   //转为double
                var H = Vector<double>.Build.DenseOfArray(R0XV0Array);
                double h = H.Norm(2);
                double p = Math.Pow(h, 2) / mu;
                double aa = -(V0.DotProduct(V0) * r0 - 2 * mu) / (mu * r0);  //半长轴a的倒数,避免奇异

                double T, x1, x2, x, z, C, S, t, derivtx, y;
                if (aa > 0) // 对于椭圆，我们用二分法
                {
                    T = 2 * Math.PI * Math.Sqrt(Math.Pow(aa, -3) / mu);
                    while (t1 - t0 > 0)
                    {
                        t1 = t1 - T;
                    }
                    t1 = t1 + T;  // t1 一定是 t0 之后一个周期内的时刻
                    x1 = 0;
                    x2 = Math.Sqrt(1 / aa) * 2 * Math.PI;
                    x = (x1 + x2) / 2;
                    z = Math.Pow(x, 2) / (1 / aa);
                    if (Math.Abs(z) < Math.Pow(10, -6))
                    {
                        C = 0.5;
                        S = 1 / 6;
                    }
                    else
                    {
                        if (z > 0)
                        {
                            C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                            S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                        }
                        else
                        {
                            C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                            S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                        }
                    }
                    while (x2 - x1 >= 0.00001)
                    {
                        x = (x1 + x2) / 2;
                        z = Math.Pow(x, 2) / (1 / aa);
                        if (Math.Abs(z) < Math.Pow(10, -6))
                        {
                            C = 0.5;
                            S = 1 / 6;
                        }
                        else
                        {
                            if (z > 0)
                            {
                                C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                                S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                            }
                            else
                            {
                                C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                                S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                            }
                        }
                        t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + R0.DotProduct(V0) / Math.Sqrt(mu) * Math.Pow(x, 2) * C + r0 * x * (1 - z * S)) + t0;
                        if (t < t1)
                        {
                            x1 = x;
                        }
                        else
                        {
                            x2 = x;
                        }
                    }
                }
                else // 对于非椭圆，我们用牛顿下山法
                {
                    derivtx = 1; // derivtx 预赋初值
                    y = Math.Sqrt(-1 / aa) * Math.Log((-2 * mu * (t1 - t0)) / (1 / aa * (R0.DotProduct(V0) + Math.Sqrt(-mu / aa) * (1 - r0 * aa)))); // x0 先用 y 表示
                    t = t1 + 5; x = y; z = Math.Pow(x, 2) / (1 / aa);
                    double xa = x;
                    if (Math.Abs(z) < Math.Pow(10, -6))
                    {
                        C = 0.5; S = 1 / 6;
                    }
                    else
                    {
                        if (z > 0)
                        {
                            C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                            S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                        }
                        else
                        {
                            C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                            S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                        }
                    }
                    t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + R0.DotProduct(V0) / Math.Sqrt(mu) * Math.Pow(x, 2) * C + r0 * x * (1 - z * S)) + t0;
                    double ta = t;
                    while (Math.Abs(t - t1) > 0.00001)
                    {
                        x = y; z = Math.Pow(x, 2) / (1 / aa);
                        if (Math.Abs(z) < Math.Pow(10, -6))
                        {
                            C = 0.5; S = 1 / 6;
                        }
                        else
                        {
                            if (z > 0)
                            {
                                C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                                S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                            }
                            else
                            {
                                C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                                S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                            }
                        }
                        t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + R0.DotProduct(V0) / Math.Sqrt(mu) * Math.Pow(x, 2) * C + r0 * x * (1 - z * S)) + t0;
                        if (Math.Abs(t - t1) > Math.Abs(ta - t1)) // 下山 
                        {
                            while (Math.Abs(t - t1) > Math.Abs(ta - t1))
                            {
                                derivtx = 2 * derivtx; y = xa - (ta - t1) / derivtx;
                                x = y;
                                z = Math.Pow(x, 2) / (1 / aa);
                                if (Math.Abs(z) < Math.Pow(10, -6))
                                {
                                    C = 0.5; S = 1 / 6;
                                }
                                else
                                {
                                    if (z > 0)
                                    {
                                        C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                                        S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                                    }
                                    else
                                    {
                                        C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                                        S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                                    }
                                }
                                t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + R0.DotProduct(V0) / Math.Sqrt(mu) * Math.Pow(x, 2) * C + r0 * x * (1 - z * S)) + t0;
                            }
                        }
                        else
                        {
                            ta = t; derivtx = 1 / Math.Sqrt(mu) * (Math.Pow(x, 2) * C + R0.DotProduct(V0) / Math.Sqrt(mu) * x * (1 - z * S) + r0 * (1 - z * C)); // dt/dx 
                            if (Math.Abs(derivtx) < 0.001)
                            {
                                derivtx = derivtx + 1;
                            }
                            y = x - (t - t1) / derivtx; xa = x;
                        }
                    }
                }

                // 求普适变量表示的状态转移函数
                double f = 1 - Math.Pow(x, 2) / r0 * C;
                double g = t1 - t0 - 1 / Math.Sqrt(mu) * Math.Pow(x, 3) * S;
                R1 = f * R0 + g * V0;
                double r1 = R1.Norm(2);
                double fv = -Math.Sqrt(mu) / (r0 * r1) * x * (1 - z * S);
                //状态转移矩阵中的 f'，r1 是 R1 的模长 
                double gv = 1 - Math.Pow(x, 2) / r1 * C;

                V1 = fv * R0 + gv * V0;
            }
        }

        public static void SolveLambertMultiRevolution(Vector<double> R1, Vector<double> R2, double t1, double t2, Vector<double> V1b, Vector<double> V2a, out double flag, out Vector<double> V1, out Vector<double> V2)
        {
            double mu = 398600.5 * Math.Pow(10, 9);
            flag = 0;
            V1 = CreateVector.Dense<double>(3, 1e10);
            V2 = CreateVector.Dense<double>(3, 1e10);

            if (Math.Abs(t1 - t2) < 0.001)
            {
                return;
            }

            double r1 = R1.L2Norm();
            double r2 = R2.L2Norm();
            double theta_min = Math.Acos(Math.Max(Math.Min(R1.DotProduct(R2) / (r1 * r2), 1), -1));
            double theta_max = 2 * Math.PI - theta_min;
            double theta = (VectorCrossProduct(R1, R2)[2] >= 0) ? theta_min : theta_max;

            if (theta <= 0.001 || theta >= 2 * Math.PI - 0.001)
            {
                return;
            }

            double A = Math.Sqrt(r1 * r2) * Math.Sin(theta) / Math.Sqrt(1 - Math.Cos(theta));
            Vector<double> VB1 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> VB2 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> VA1 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> VA2 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> V11 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> V12 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> V21 = CreateVector.Dense<double>(3, 1e10);
            Vector<double> V22 = CreateVector.Dense<double>(3, 1e10);

            double z1 = Math.Pow(2 * Math.PI, 2);
            double z2 = Math.Pow(4 * Math.PI, 2);
            double epsino = 1;
            double k = 0.618;
            double L = z1 + (1 - k) * (z2 - z1);
            double R = z1 + k * (z2 - z1);

            while (epsino >= 0.000001)
            {
                double CL = (1 - Math.Cos(Math.Sqrt(L))) / L;
                double SL = (Math.Sqrt(L) - Math.Sin(Math.Sqrt(L))) / Math.Pow(Math.Sqrt(L), 3);
                double CR = (1 - Math.Cos(Math.Sqrt(R))) / R;
                double SR = (Math.Sqrt(R) - Math.Sin(Math.Sqrt(R))) / Math.Pow(Math.Sqrt(R), 3);
                double yL = r1 + r2 - A * (1 - L * SL) / Math.Sqrt(CL);
                double yR = r1 + r2 - A * (1 - R * SR) / Math.Sqrt(CR);
                double xL = Math.Sqrt(yL / CL);
                double xR = Math.Sqrt(yR / CR);
                double tL = 1 / Math.Sqrt(mu) * (Math.Pow(xL, 3) * SL + A * Math.Sqrt(yL)) + t1;
                double tR = 1 / Math.Sqrt(mu) * (Math.Pow(xR, 3) * SR + A * Math.Sqrt(yR)) + t1;

                if (tL > tR)
                {
                    epsino = R - L;
                    z1 = L;
                    L = R;
                    R = z1 + k * (z2 - z1);
                }
                else
                {
                    epsino = R - L;
                    z2 = R;
                    R = L;
                    L = z1 + (1 - k) * (z2 - z1);
                }
            }

            double zBmin = (L + R) / 2;
            double tBmin = t1;
            double y1m = 1.0, z1m = 1.0, C1m = 1.0, S1m = 1.0, x1m = 0, t1m = 0;

            if (t2 > tBmin)
            {
                // 多圈判断 case1
                double z1L = Math.Pow(2 * Math.PI, 2);
                double z1R = zBmin;

                while (z1R - z1L >= 0.0000001)
                {
                    z1m = (z1L + z1R) / 2;
                    C1m = (1 - Math.Cos(Math.Sqrt(z1m))) / z1m;
                    S1m = (Math.Sqrt(z1m) - Math.Sin(Math.Sqrt(z1m))) / Math.Pow(Math.Sqrt(z1m), 3);
                    y1m = r1 + r2 - A * (1 - z1m * S1m) / Math.Sqrt(C1m);
                    x1m = Math.Sqrt(y1m / C1m);
                    t1m = 1 / Math.Sqrt(mu) * (Math.Pow(x1m, 3) * S1m + A * Math.Sqrt(y1m)) + t1;

                    if (t1m > t2)
                    {
                        z1L = z1m;
                    }
                    else
                    {
                        z1R = z1m;
                    }
                }

                double z2L = zBmin, z2R = Math.Pow(4 * Math.PI, 2);
                double z2m = 1.0, C2m = 1.0, S2m = 1.0, y2m = 1.0, x2m = 0, t2m = 0;

                while (z2R - z2L >= 0.00001)
                {
                    z2m = (z2L + z2R) / 2;
                    C2m = (1 - Math.Cos(Math.Sqrt(z2m))) / z2m;
                    S2m = (Math.Sqrt(z2m) - Math.Sin(Math.Sqrt(z2m))) / Math.Pow(Math.Sqrt(z2m), 3);
                    y2m = r1 + r2 - A * (1 - z2m * S2m) / Math.Sqrt(C2m);
                    x2m = Math.Sqrt(y2m / C2m);
                    t2m = 1 / Math.Sqrt(mu) * (Math.Pow(x2m, 3) * S2m + A * Math.Sqrt(y2m)) + t1;

                    if (t2m > t2)
                    {
                        z2R = z2m;
                    }
                    else
                    {
                        z2L = z2m;
                    }
                }

                if (Math.Abs(theta - Math.PI) < 0.001)
                {
                    Vector<double> nx = R1 / r1;
                    Vector<double> ny = VectorCrossProduct(VectorCrossProduct(R1, V1b), R1).Normalize(2);
                    V11 = -Math.Sqrt(mu) * (1 - z1m * S1m) / Math.Sqrt(C1m * (r1 + r2)) * nx + Math.Sqrt(2 * mu * r2) / Math.Sqrt(r1 * (r1 + r2)) * ny;
                    V21 = -Math.Sqrt(mu) * (1 - z1m * S1m) / Math.Sqrt(C1m * (r1 + r2)) * nx - Math.Sqrt(2 * mu * r1) / Math.Sqrt(r2 * (r1 + r2)) * ny;
                }
                else
                {
                    double f = 1 - y1m / r1;
                    double g = A * Math.Sqrt(y1m / mu);
                    double gv = 1 - y1m / r2;
                    VB1 = (R2 - f * R1) / g;
                    VB2 = (gv * R2 - R1) / g;
                }

                if (Math.Abs(theta - Math.PI) < 0.001)
                {
                    Vector<double> nx = R1 / r1;
                    Vector<double> ny = VectorCrossProduct(VectorCrossProduct(R1, V1b), R1).Normalize(2);
                    V12 = -Math.Sqrt(mu) * (1 - z2m * S2m) / Math.Sqrt(C2m * (r1 + r2)) * nx + Math.Sqrt(2 * mu * r2) / Math.Sqrt(r1 * (r1 + r2)) * ny;
                    V22 = -Math.Sqrt(mu) * (1 - z2m * S2m) / Math.Sqrt(C2m * (r1 + r2)) * nx - Math.Sqrt(2 * mu * r1) / Math.Sqrt(r2 * (r1 + r2)) * ny;
                }
                else
                {
                    double f = 1 - y2m / r1;
                    double g = A * Math.Sqrt(y2m / mu);
                    double gv = 1 - y2m / r2;
                    VA1 = (R2 - f * R1) / g;
                    VA2 = (gv * R2 - R1) / g;
                }

                if ((V11 - V1b).L2Norm() + (V21 - V2a).L2Norm() < (V12 - V1b).L2Norm() + (V22 - V2a).L2Norm())
                {
                    VB1 = V11; // case1 的 t1 脉冲
                    VB2 = V21; // case1 的 t2 脉冲
                }
                else
                {
                    VB1 = V12;
                    VB2 = V22;
                }
            }

            if (t2 == tBmin)
            {
                double CB = (1 - Math.Cos(Math.Sqrt(zBmin))) / zBmin;
                double SB = (Math.Sqrt(zBmin) - Math.Sin(Math.Sqrt(zBmin))) / Math.Pow(Math.Sqrt(zBmin), 3);
                double yB = r1 + r2 - A * (1 - zBmin * SB) / Math.Sqrt(CB);
                double xB = Math.Sqrt(yB / CB);
                double tB = 1 / Math.Sqrt(mu) * (Math.Pow(xB, 3) * SB + A * Math.Sqrt(yB)) + t1;
                // 多圈判断 case2
                if (Math.Abs(theta - Math.PI) < 0.001)
                {
                    Vector<double> nx = R1 / r1;
                    Vector<double> ny = VectorCrossProduct(VectorCrossProduct(R1, V1b), R1).Normalize(2);
                    VB1 = -Math.Sqrt(mu) * (1 - zBmin * SB) / Math.Sqrt(CB * (r1 + r2)) * nx + Math.Sqrt(2 * mu * r2) / Math.Sqrt(r1 * (r1 + r2)) * ny;
                    VB2 = -Math.Sqrt(mu) * (1 - zBmin * SB) / Math.Sqrt(CB * (r1 + r2)) * nx - Math.Sqrt(2 * mu * r1) / Math.Sqrt(r2 * (r1 + r2)) * ny;
                }
                else
                {
                    double f = 1 - yB / r1;
                    double g = A * Math.Sqrt(yB / mu);
                    double gv = 1 - yB / r2;
                    VB1 = (R2 - f * R1) / g;
                    VB2 = (gv * R2 - R1) / g;
                }
            }
            // %case2结束
            // %多圈情况脉冲求解结束


            // 单圈判断
            double zz = 0;    //z的初值
            double za = zz + 10;
            double z = zz;
            double C = 0.5;
            double S = 1.0 / 6.0;
            double y = r1 + r2 - A * (1 - zz * S) / Math.Sqrt(C);
            double x = Math.Sqrt(y / C);
            double t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + A * Math.Sqrt(y)) + t1;

            if (t < t2)    // 对z在0-4pi^2用二分法
            {
                double zmin = 0;
                double zmax = Math.Pow(2 * Math.PI, 2);

                while (Math.Abs(t - t2) >= 1e-9)
                {
                    z = (zmax + zmin) / 2;
                    if (Math.Abs(z - zmin) < 1e-9)
                    {
                        break;
                    }
                    C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                    S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                    y = r1 + r2 - A * (1 - z * S) / Math.Sqrt(C);
                    x = Math.Sqrt(y / C);
                    t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + A * Math.Sqrt(y)) + t1;
                    if (t > t2)
                    {
                        zmax = z;
                    }
                    else
                    {
                        zmin = z;
                    }
                }
            }
            else
            {
                int times = 0;
                double derivtz = 1; // derivtz 预赋初值
                while (Math.Abs(t - t2) >= 1e-9)
                {
                    z = zz;
                    if (Math.Abs(z - za) < 1e-9)
                    {
                        break;
                    }
                    if (Math.Abs(z) < 0.0001)
                    {
                        C = 0.5;
                        S = 1.0 / 6.0;
                    }
                    else
                    {
                        if (z > 0)
                        {
                            C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                            S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                        }
                        else
                        {
                            C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                            S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                        }
                    }
                    y = r1 + r2 - A * (1 - z * S) / Math.Sqrt(C);
                    // 对y为负时处理
                    int times2 = 0;
                    Random rand = new Random();
                    while (y < 0)
                    {
                        double randk = 2 * rand.NextDouble();
                        z = za - (t - t2) / (randk * derivtz);
                        derivtz = randk * derivtz;
                        if (Math.Abs(z) < 0.0001)
                        {
                            C = 0.5;
                            S = 1.0 / 6.0;
                        }
                        else
                        {
                            if (z > 0)
                            {
                                C = (1 - Math.Cos(Math.Sqrt(z))) / z;
                                S = (Math.Sqrt(z) - Math.Sin(Math.Sqrt(z))) / Math.Pow(Math.Sqrt(z), 3);
                            }
                            else
                            {
                                C = (Math.Cosh(Math.Sqrt(-z)) - 1) / (-z);
                                S = (Math.Sinh(Math.Sqrt(-z)) - Math.Sqrt(-z)) / Math.Pow(Math.Sqrt(-z), 3);
                            }
                        }
                        y = r1 + r2 - A * (1 - z * S) / Math.Sqrt(C);
                        za = z;
                        times2++;
                        if (times2 > 100000)
                        {
                            times = 100000 + 1;
                            break;
                        }
                    }

                    x = Math.Sqrt(y / C);
                    t = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * S + A * Math.Sqrt(y)) + t1;
                    // 求导
                    double DSz = 1 / (2 * z) * (C - 3 * S);
                    double DCz = 1 / (2 * z) * (1 - z * S - 2 * C);
                    if (Math.Abs(z) < 0.0001)
                    {
                        derivtz = 1 / Math.Sqrt(mu) * (Math.Sqrt(2) / 40 * Math.Pow(y, 1.5) + A / 8 * (Math.Sqrt(y) + A * Math.Sqrt(1 / (2 * y))));
                    }
                    else
                    {
                        derivtz = 1 / Math.Sqrt(mu) * (Math.Pow(x, 3) * (DSz - 3 * S * DCz / (2 * C)) + A / 8 * (3 * S * Math.Sqrt(y) / C + A / x));
                    }
                    // 下一步迭代的z先存在zz
                    zz = z - (t - t2) / derivtz;
                    if (zz > Math.Pow(2 * Math.PI, 2))
                    {
                        zz = 7 * rand.NextDouble() + 32;
                    }
                    za = z; // 存储这一步的 z
                    times++;
                    if (times > 1e5)
                    {
                        break;
                    }
                }
            }
            // 求两端点速度
            // 对于真近点角差为pi定义非奇异公式
            if (Math.Abs(theta - Math.PI) < 0.001)
            {
                Vector<double> nx = R1 / r1;
                // Vector<double> nyvector = VectorCrossProduct(VectorCrossProduct(R1, V1b), R1);
                // double nynorm = nyvector.L2Norm();
                // Vector<double> ny = vectorToNormalize / nynorm;
                Vector<double> ny = VectorCrossProduct(VectorCrossProduct(R1, V1b), R1).Normalize(2);
                VA1 = -Math.Sqrt(mu) * (1 - z * S) / Math.Sqrt(C * (r1 + r2)) * nx + Math.Sqrt(2 * mu * r2) / Math.Sqrt(r1 * (r1 + r2)) * ny;
                VA2 = -Math.Sqrt(mu) * (1 - z * S) / Math.Sqrt(C * (r1 + r2)) * nx - Math.Sqrt(2 * mu * r1) / Math.Sqrt(r2 * (r1 + r2)) * ny;
            }
            else
            {
                double f = 1 - y / r1;
                double g = A * Math.Sqrt(y / mu);
                double gv = 1 - y / r2;
                VA1 = (R2 - f * R1) / g;
                VA2 = (gv * R2 - R1) / g;
            }
            // 单圈求解结束
            // 比较单圈和多圈结果
            if ((VA1 - V1b).L2Norm() + (VA2 - V2a).L2Norm() <= (VB1 - V1b).L2Norm() + (VB2 - V2a).L2Norm())
            {
                V1 = VA1;
                V2 = VA2;
                flag = 1;
            }
            else
            {
                V1 = VB1;
                V2 = VB2;
                flag = 2;
            }
        }

        public static Vector<double> VectorCrossProduct(Vector<double> vector1, Vector<double> vector2)
        {
            return CreateVector.DenseOfArray(new double[]
            {
                vector1[1] * vector2[2] - vector1[2] * vector2[1],
                vector1[2] * vector2[0] - vector1[0] * vector2[2],
                vector1[0] * vector2[1] - vector1[1] * vector2[0]
            });
        }

        public static void TriImpRendez(double tA, double tB, double tC, double deltav1,
                                       double alpha, double beta, Vector<double> RT0, Vector<double> VT0,
                                       Vector<double> RP0, Vector<double> VP0, double t0, out double J)
        {
            // Sort the time variables to get t1, t2 and tf
            double[] tsq = { tA, tB, tC };
            Array.Sort(tsq);
            double t1 = tsq[0];
            double t2 = tsq[1];
            double tf = tsq[2];
            // 定义常数
            double re = 6371.11 * Math.Pow(10, 3);
            double hinf = re + 150 * Math.Pow(10, 3);
            double mu = 398600.5 * Math.Pow(10, 9);
            // 评价函数
            J = 1e10;

            // 求解最终状态
            CombineSolvePropagation(tf, RT0, VT0, t0, out Vector<double> RTf, out Vector<double> VTf);
            CombineSolvePropagation(t1, RP0, VP0, t0, out Vector<double> RP1, out Vector<double> VP1b);

            // Compute the DeltaV1 vector
            Vector<double> DeltaV1 = Vector<double>.Build.DenseOfArray(new double[] {
                deltav1 * Math.Cos(alpha) * Math.Cos(beta),
                deltav1 * Math.Cos(alpha) * Math.Sin(beta),
                deltav1 * Math.Sin(alpha)
            });

            Vector<double> VP1a = VP1b + DeltaV1;

            CombineSolvePropagation(t2, RP1, VP1a, t1, out Vector<double> RP2, out Vector<double> VP2b);


            // 判断高度约束是否满足
            double rp1 = RP1.L2Norm();
            double vp1a = VP1a.L2Norm();
            double rp2 = RP2.L2Norm();
            double vp2b = VP2b.L2Norm();
            double aa = -(vp1a * vp1a * rp1 - 2 * mu) / (mu * rp1);

            Vector<double> vector_e = (1 / mu) * ((vp1a * vp1a - mu / rp1) * RP1 - RP1.DotProduct(VP1a) * VP1a);
            double e = vector_e.L2Norm();

            double v_cos_1 = vector_e.DotProduct(RP1) / (e * rp1);
            v_cos_1 = (v_cos_1 > 1) ? 1 : v_cos_1;
            v_cos_1 = (v_cos_1 < -1) ? -1 : v_cos_1;

            double v_cos_2 = vector_e.DotProduct(RP2) / (e * rp2);
            v_cos_2 = (v_cos_2 > 1) ? 1 : v_cos_2;
            v_cos_2 = (v_cos_2 < -1) ? -1 : v_cos_2;

            double hmin1;
            if (aa > 0) //对于该段轨道是椭圆弧的情况
            {
                double T = 2 * Math.PI * Math.Sqrt(1 / (aa * aa * aa * mu));
                if (t2 - t1 >= T) //飞行时间大于一个周期，近地点为最低点
                {
                    hmin1 = 1 / aa * (1 - e);
                }
                else //飞行时间小于一个周期
                {   //真近点角在0到pi还是在pi到2pi的范围内
                    double theta1 = (RP1.DotProduct(VP1a) >= 0) ? Math.Acos(v_cos_1) : 2 * Math.PI - Math.Acos(v_cos_1);
                    //真近点角在0到pi还是在pi到2pi的范围内
                    double theta2 = (RP2.DotProduct(VP2b) >= 0) ? Math.Acos(v_cos_2) : 2 * Math.PI - Math.Acos(v_cos_2);
                    // 是否经过近地点
                    hmin1 = (theta1 <= theta2) ? Math.Min(rp1, rp2) : 1 / aa * (1 - e);
                }
            }
            else  //对于轨道非椭圆的情况
            {
                //                              真近点角在0到某个正数之间       真近点角在某个负数到0之间  这句先不改，回头看看需不要加2pi
                double theta1 = (RP1.DotProduct(VP1a) >= 0) ? Math.Acos(v_cos_1) : -Math.Acos(v_cos_1);
                double theta2 = (RP2.DotProduct(VP2b) >= 0) ? Math.Acos(v_cos_2) : -Math.Acos(v_cos_2);
                double dtheta = theta2 - theta1;
                double h = Math.Pow(vector_e.L2Norm(), 2) / (2 * mu) * (1 / Math.Cos(dtheta / 2) - 1);
                if (theta2 <= 0)
                {
                    hmin1 = rp2;
                }
                else if (theta1 >= 0)
                {
                    hmin1 = rp1;
                }
                else
                {
                    hmin1 = 1 / aa * (1 - e);
                }
            }
            if (hmin1 < hinf)//判断轨道高度约束结束
            {
                return;
            }


            SolveLambertMultiRevolution(RP2, RTf, t2, tf, VP2b, VTf, out double flag, out Vector<double> VP2a, out Vector<double> VPfb);
            if (flag == 0)
            {
                return;
            }
            rp2 = RP2.L2Norm();
            double vp2a = VP2a.L2Norm();
            double rtf = RTf.L2Norm();
            double vpfb = VPfb.L2Norm();
            aa = -(vp2a * vp2a * rp2 - 2 * mu) / (mu * rp2);
            //偏心率矢量
            vector_e = 1 / mu * ((vp2a * vp2a - mu / rp2) * RP2 - RP2.DotProduct(VP2a) * VP2a);
            e = vector_e.L2Norm();

            v_cos_1 = vector_e.DotProduct(RP2) / (e * rp2);
            v_cos_1 = (v_cos_1 > 1) ? 1 : v_cos_1;
            v_cos_1 = (v_cos_1 < -1) ? -1 : v_cos_1;

            v_cos_2 = vector_e.DotProduct(RTf) / (e * rtf);
            v_cos_2 = (v_cos_2 > 1) ? 1 : v_cos_2;
            v_cos_2 = (v_cos_2 < -1) ? -1 : v_cos_2;
            double hmin2;
            if (aa > 0) //对于该段轨道是椭圆弧的情况
            {
                double T = 2 * Math.PI * Math.Sqrt(1 / (aa * aa * aa * mu));
                if (tf - t2 >= T) //飞行时间大于一个周期，近地点为最低点
                {
                    hmin2 = 1 / aa * (1 - e);
                }
                else  //飞行时间小于一个周期
                {
                    //真近点角在0到pi还是在pi到2pi的范围内
                    double theta1 = (RP2.DotProduct(VP2a) >= 0) ? Math.Acos(v_cos_1) : 2 * Math.PI - Math.Acos(v_cos_1);
                    double theta2 = (RTf.DotProduct(VPfb) >= 0) ? Math.Acos(v_cos_2) : 2 * Math.PI - Math.Acos(v_cos_2);
                    // 有没有经过近地点
                    hmin2 = (theta1 <= theta2) ? Math.Min(rp2, rtf) : 1 / aa * (1 - e);
                }
            }
            else //对于轨道非椭圆的情况
            {
                double theta1 = (RP2.DotProduct(VP2a) >= 0) ? Math.Acos(v_cos_1) : -Math.Acos(v_cos_1);
                double theta2 = (RTf.DotProduct(VPfb) >= 0) ? Math.Acos(v_cos_2) : -Math.Acos(v_cos_2);
                if (theta2 <= 0)
                {
                    hmin2 = rtf;
                }
                else if (theta1 >= 0)
                {
                    hmin2 = rp2;
                }
                else
                {
                    hmin2 = 1 / aa * (1 - e);
                }
            }
            if (hmin2 < hinf)
            {
                return;
            }
            Vector<double> DeltaV2 = VP2a - VP2b;
            double deltav2 = DeltaV2.L2Norm();
            Vector<double> DeltaV3 = VTf - VPfb;
            double deltav3 = DeltaV3.L2Norm();
            J = deltav1 + deltav2 + deltav3;

            if (double.IsNaN(J) || double.IsInfinity(J))
            {
                J = 1e10;
            }

        }

        public static void SolveImp(double tA, double tB, double tC, double deltav1,
                               double alpha, double beta, Vector<double> RT0, Vector<double> VT0,
                               Vector<double> RP0, Vector<double> VP0, double t0, out Vector<double> DeltaV1, out Vector<double> DeltaV2, out Vector<double> DeltaV3)
        {
            // Sort the time variables to get t1, t2 and tf
            double[] tsq = { tA, tB, tC };
            Array.Sort(tsq);
            double t1 = tsq[0];
            double t2 = tsq[1];
            double tf = tsq[2];
            // 定义常数
            double re = 6371.11 * Math.Pow(10, 3);
            double hinf = re + 150 * Math.Pow(10, 3);
            double mu = 398600.5 * Math.Pow(10, 9);

            // 求解最终状态
            CombineSolvePropagation(tf, RT0, VT0, t0, out Vector<double> RTf, out Vector<double> VTf);
            CombineSolvePropagation(t1, RP0, VP0, t0, out Vector<double> RP1, out Vector<double> VP1b);

            // Compute the DeltaV1 vector
            DeltaV1 = Vector<double>.Build.DenseOfArray(new double[] {
                deltav1 * Math.Cos(alpha) * Math.Cos(beta),
                deltav1 * Math.Cos(alpha) * Math.Sin(beta),
                deltav1 * Math.Sin(alpha)
            });

            Vector<double> VP1a = VP1b + DeltaV1;

            CombineSolvePropagation(t2, RP1, VP1a, t1, out Vector<double> RP2, out Vector<double> VP2b);


            // 判断高度约束是否满足
            double rp1 = RP1.L2Norm();
            double vp1a = VP1a.L2Norm();
            double rp2 = RP2.L2Norm();
            double vp2b = VP2b.L2Norm();
            double aa = -(vp1a * vp1a * rp1 - 2 * mu) / (mu * rp1);

            Vector<double> vector_e = (1 / mu) * ((vp1a * vp1a - mu / rp1) * RP1 - RP1.DotProduct(VP1a) * VP1a);
            double e = vector_e.L2Norm();

            double v_cos_1 = vector_e.DotProduct(RP1) / (e * rp1);
            v_cos_1 = (v_cos_1 > 1) ? 1 : v_cos_1;
            v_cos_1 = (v_cos_1 < -1) ? -1 : v_cos_1;

            double v_cos_2 = vector_e.DotProduct(RP2) / (e * rp2);
            v_cos_2 = (v_cos_2 > 1) ? 1 : v_cos_2;
            v_cos_2 = (v_cos_2 < -1) ? -1 : v_cos_2;

            double hmin1;
            if (aa > 0) //对于该段轨道是椭圆弧的情况
            {
                double T = 2 * Math.PI * Math.Sqrt(1 / (aa * aa * aa * mu));
                if (t2 - t1 >= T) //飞行时间大于一个周期，近地点为最低点
                {
                    hmin1 = 1 / aa * (1 - e);
                }
                else //飞行时间小于一个周期
                {   //真近点角在0到pi还是在pi到2pi的范围内
                    double theta1 = (RP1.DotProduct(VP1a) >= 0) ? Math.Acos(v_cos_1) : 2 * Math.PI - Math.Acos(v_cos_1);
                    //真近点角在0到pi还是在pi到2pi的范围内
                    double theta2 = (RP2.DotProduct(VP2b) >= 0) ? Math.Acos(v_cos_2) : 2 * Math.PI - Math.Acos(v_cos_2);
                    // 是否经过近地点
                    hmin1 = (theta1 <= theta2) ? Math.Min(rp1, rp2) : 1 / aa * (1 - e);
                }
            }
            else  //对于轨道非椭圆的情况
            {
                //                              真近点角在0到某个正数之间       真近点角在某个负数到0之间  这句先不改，回头看看需不要加2pi
                double theta1 = (RP1.DotProduct(VP1a) >= 0) ? Math.Acos(v_cos_1) : -Math.Acos(v_cos_1);
                double theta2 = (RP2.DotProduct(VP2b) >= 0) ? Math.Acos(v_cos_2) : -Math.Acos(v_cos_2);
                double dtheta = theta2 - theta1;
                double h = Math.Pow(vector_e.L2Norm(), 2) / (2 * mu) * (1 / Math.Cos(dtheta / 2) - 1);
                if (theta2 <= 0)
                {
                    hmin1 = rp2;
                }
                else if (theta1 >= 0)
                {
                    hmin1 = rp1;
                }
                else
                {
                    hmin1 = 1 / aa * (1 - e);
                }
            }
            //if (hmin1 < hinf)//判断轨道高度约束结束
            //{
            //    return;
            //}


            SolveLambertMultiRevolution(RP2, RTf, t2, tf, VP2b, VTf, out double flag, out Vector<double> VP2a, out Vector<double> VPfb);
            //if (flag == 0)
            //{
            //    return;
            //}
            rp2 = RP2.L2Norm();
            double vp2a = VP2a.L2Norm();
            double rtf = RTf.L2Norm();
            double vpfb = VPfb.L2Norm();
            aa = -(vp2a * vp2a * rp2 - 2 * mu) / (mu * rp2);
            //偏心率矢量
            vector_e = 1 / mu * ((vp2a * vp2a - mu / rp2) * RP2 - RP2.DotProduct(VP2a) * VP2a);
            e = vector_e.L2Norm();

            v_cos_1 = vector_e.DotProduct(RP2) / (e * rp2);
            v_cos_1 = (v_cos_1 > 1) ? 1 : v_cos_1;
            v_cos_1 = (v_cos_1 < -1) ? -1 : v_cos_1;

            v_cos_2 = vector_e.DotProduct(RTf) / (e * rtf);
            v_cos_2 = (v_cos_2 > 1) ? 1 : v_cos_2;
            v_cos_2 = (v_cos_2 < -1) ? -1 : v_cos_2;
            double hmin2;
            if (aa > 0) //对于该段轨道是椭圆弧的情况
            {
                double T = 2 * Math.PI * Math.Sqrt(1 / (aa * aa * aa * mu));
                if (tf - t2 >= T) //飞行时间大于一个周期，近地点为最低点
                {
                    hmin2 = 1 / aa * (1 - e);
                }
                else  //飞行时间小于一个周期
                {
                    //真近点角在0到pi还是在pi到2pi的范围内
                    double theta1 = (RP2.DotProduct(VP2a) >= 0) ? Math.Acos(v_cos_1) : 2 * Math.PI - Math.Acos(v_cos_1);
                    double theta2 = (RTf.DotProduct(VPfb) >= 0) ? Math.Acos(v_cos_2) : 2 * Math.PI - Math.Acos(v_cos_2);
                    // 有没有经过近地点
                    hmin2 = (theta1 <= theta2) ? Math.Min(rp2, rtf) : 1 / aa * (1 - e);
                }
            }
            else //对于轨道非椭圆的情况
            {
                double theta1 = (RP2.DotProduct(VP2a) >= 0) ? Math.Acos(v_cos_1) : -Math.Acos(v_cos_1);
                double theta2 = (RTf.DotProduct(VPfb) >= 0) ? Math.Acos(v_cos_2) : -Math.Acos(v_cos_2);
                if (theta2 <= 0)
                {
                    hmin2 = rtf;
                }
                else if (theta1 >= 0)
                {
                    hmin2 = rp2;
                }
                else
                {
                    hmin2 = 1 / aa * (1 - e);
                }
            }
            //if (hmin2 < hinf)
            //{
            //    return;
            //}
            DeltaV2 = VP2a - VP2b;
            DeltaV3 = VTf - VPfb;

        }

        public static void Randperm(int D, out int[] ar)
        {
            Random randperm = new Random(); // 创建一个随机数生成器
            ar = Enumerable.Range(1, D).ToArray(); // 创建包含 1 到 D 的有序数组
                                                   // 使用 Fisher–Yates 洗牌算法打乱数组
            for (int i = 0; i < D - 1; i++)
            {
                int j = randperm.Next(i, D);
                int temp = ar[i];
                ar[i] = ar[j];
                ar[j] = temp;
            }
        }

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

        // public static void CLPSO_trp(int Dimension, Vector<double> VRmin, Vector<double> VRmax, int Max_Gen, int Particle_Number, Vector<double> RT0, Vector<double> VT0, Vector<double> RP0, Vector<double> VP0, double t0, out Vector<double> gbest, out double gbestval, out double[] gb_record, out double[] ii)
        public static void PSO_trp(Func<Vector<double>, double> fhd, int Dimension, Vector<double> Rmin, Vector<double> Rmax, int Max_Gen, int Particle_Number, out Vector<double> gbest, out double gbestval)
        {
            int ps = Particle_Number; // 初始种群个数
            int D = Dimension;        // 空间维数
            int me = Max_Gen;         // 最大迭代次数

            if (Rmin.Count == 1)
            {
                Rmin = Vector<double>.Build.Dense(D, Rmin[0]);
                Rmax = Vector<double>.Build.Dense(D, Rmax[0]);
            }

            Vector<double> mv = 0.2 * (Rmax - Rmin);
            Matrix<double> Rmin_M = Matrix<double>.Build.Dense(ps, D, (i, j) => Rmin[j]);
            Matrix<double> Rmax_M = Matrix<double>.Build.Dense(ps, D, (i, j) => Rmax[j]);
            Matrix<double> Vmin_M = -mv.ToColumnMatrix();
            Matrix<double> Vmax_M = -Vmin_M;

            Vector<double> w = Vector<double>.Build.Dense(me, i => 0.9 - ((double)i * 0.7 / (me - 1)));
            double c1 = 0.8; // 自我学习因子
            double c2 = 1.49; // 群体学习因子

            Random rand = new Random();
            Matrix<double> pos = Matrix<double>.Build.Dense(ps, D, (i, j) => Rmin[j] + rand.NextDouble() * (Rmax[j] - Rmin[j]));
            Matrix<double> vel = Matrix<double>.Build.Dense(ps, D, (i, j) => Vmin_M[i, j] + rand.NextDouble() * (Vmax_M[i, j] - Vmin_M[i, j]));

            List<Vector<double>> positions = pos.EnumerateRows().Select(row => row.Clone()).ToList();
            List<Vector<double>> velocities = vel.EnumerateRows().Select(row => row.Clone()).ToList();
            List<double> fitnesses = positions.Select(fhd).ToList();
            List<Vector<double>> pBestPositions = positions.Select(row => row.Clone()).ToList();
            List<double> pBestFitnesses = fitnesses.ToList();

            gbestval = pBestFitnesses.Min();
            int minIndex = pBestFitnesses.IndexOf(gbestval);
            gbest = pBestPositions[minIndex].Clone();
            Vector<double> gbestLocal = gbest.Clone(); // 使用本地变量存储 gbest

            for (int i = 1; i < me; ++i)
            {
                for (int k = 0; k < ps; ++k)
                {
                    velocities[k] = w[i] * velocities[k] +
                                   c1 * Vector<double>.Build.Dense(D, j => (pBestPositions[k][j] - positions[k][j]) * rand.NextDouble()) +
                                   c2 * Vector<double>.Build.Dense(D, j => (gbestLocal[j] - positions[k][j]) * rand.NextDouble());

                    velocities[k] = Vector<double>.Build.Dense(D, j => Math.Max(-mv[j], Math.Min(mv[j], velocities[k][j])));
                    positions[k] = positions[k] + velocities[k];

                    positions[k] = Vector<double>.Build.Dense(D, j => Math.Max(Rmin_M[k, j], Math.Min(Rmax_M[k, j], positions[k][j])));

                    fitnesses[k] = fhd(positions[k]);

                    if (fitnesses[k] <= pBestFitnesses[k])
                    {
                        pBestPositions[k] = positions[k].Clone();
                        pBestFitnesses[k] = fitnesses[k];
                    }

                    if (pBestFitnesses[k] < gbestval)
                    {
                        gbest = pBestPositions[k].Clone();
                        gbestval = pBestFitnesses[k];
                    }
                }
            }
        }


        //////////测试 TriImpRendez
        //private static void Main(string[] args)
        //{
        //    // 目标函数，维度，变量下界，变量上界，最大代数，最大计算目标函数次数，最大迟滞代数，精度      
        //    double mu = 398600.5E9;


        //    // space shuttle orbital elemet
        //    double aP = 42.166E6;
        //    double eP = 0;
        //    double iP = 0 * Math.PI / 180;
        //    double RAANP = 0 * Math.PI / 180;
        //    double wP = 0 * Math.PI / 180;
        //    double thetaP = -0.0118578949864820;
        //    double pP = aP * (1 - eP * eP);
        //    double TP = 2 * Math.PI * Math.Sqrt(aP * aP * aP / mu);

        //    // ISS  orbital elemet
        //    double aT = 42.166E6;
        //    double eT = 0;
        //    double iT = 0 * Math.PI / 180;
        //    double RAANT = 0 * Math.PI / 180;
        //    double wT = 0 * Math.PI / 180;
        //    double thetaT = -0.0011857894986482;
        //    double pT = aT * (1 - eT * eT);
        //    double TT = 2 * Math.PI * Math.Sqrt(aT * aT * aT / mu);

        //    // space shuttle orbital elemet in to cartesian R0 V0
        //    ElementsToRV(aP, eP, iP, RAANP, wP, thetaP, out Vector<double> RP0, out Vector<double> VP0);
        //    // ISS orbital elemetin to cartesian R0 V0
        //    ElementsToRV(aT, eT, iT, RAANT, wT, thetaT, out Vector<double> RT0, out Vector<double> VT0);
        //    double tA, tB, tC, deltav1, alpha, beta, t0;
        //    tA = 999.999744559997;
        //    tB = 0.000146089963668849;
        //    tC = 0.000250651898709001;
        //    deltav1 = 10.7568974605388;
        //    alpha = 0.000225567332278236;
        //    beta = 1.63222204069228;
        //    t0 = 0;
        //    TriImpRendez(tA, tB, tC, deltav1, alpha, beta, RT0, VT0, RP0, VP0, t0, out double J1);
        //    //////下一步是写解三脉冲
        //    Console.WriteLine($"J1 = {J1}");
        //    Console.ReadKey();
        //}

        public class Particle
        {
            public Vector<double> Position;
            public Vector<double> Velocity;
            public Vector<double> PBest;
            public double Fitness;
            public double PBestFitness;
        }

        private static void Main(string[] args)
        {
            
            double mu = 398600.5E9;
            //double re = 6371.11E3;

            Vector<double> T_oe = Vector<double>.Build.Dense(new double[] { 42.166e6, 0.0, 1 * Math.PI / 180, 0 * Math.PI / 180, 0 * Math.PI / 180, 0 * Math.PI / 180 });
            // Start distance 滞后500km
            double St_dis = -500e3;
            // Port distance 滞后10km
            double Pt_dis = -10e3;  
            // 初始时刻
            double t0 = 0;
            // 搜索上限
            double T_search = 3000;

            //-------------------------------------------- 以上为原程序输入-------------------------------------

            // 出发点航天器 orbital elemet
            double aP = T_oe[0];
            double eP = T_oe[1];
            double iP = T_oe[2];
            double RAANP = T_oe[3];
            double wP = T_oe[4];
            double thetaP = T_oe[5] + St_dis/ aP;
            double pP = aP * (1 - eP * eP);
            double TP = 2 * Math.PI * Math.Sqrt(aP * aP * aP / mu);

            // 到达点航天器 orbital elemet
            double aT = T_oe[0];
            double eT = T_oe[1];
            double iT = T_oe[2];
            double RAANT = T_oe[3];
            double wT = T_oe[4];
            double thetaT = T_oe[5] + Pt_dis / aP;
            double pT = aT * (1 - eT * eT);
            double TT = 2 * Math.PI * Math.Sqrt(aT * aT * aT / mu);

            // space shuttle orbital elemetin to cartesian R0 V0
            ElementsToRV(aP, eP, iP, RAANP, wP, thetaP, out Vector<double> RP0, out Vector<double> VP0);
            // ISS orbital elemetin to cartesian R0 V0
            ElementsToRV(aT, eT, iT, RAANT, wT, thetaT, out Vector<double> RT0, out Vector<double> VT0);



            // 维数     同时PSO代码开始于此
            int Dimension = 6;
            // 设置上下限
            Matrix<double> Rmin = DenseMatrix.OfArray(new double[,] { { 0, 0, 0, 0, -Math.PI / 2, 0 } });
            Matrix<double> Rmax = DenseMatrix.OfArray(new double[,] { { T_search, T_search, T_search, 10000, Math.PI / 2, 2 * Math.PI } });
            int Max_Gen = 1000;
            int Particle_Number = 12;

            int ps = Particle_Number;
            int D = Dimension;
            int me = Max_Gen;

            if (Rmin.RowCount == 1)
            {
                Rmin = DenseMatrix.Create(ps, D, (i, j) => Rmin[0, j]);
                Rmax = DenseMatrix.Create(ps, D, (i, j) => Rmax[0, j]);
            }

            Matrix<double> mv = 0.2 * (Rmax - Rmin);
            Matrix<double> Rmin_M = Rmin;
            Matrix<double> Rmax_M = Rmax;
            Matrix<double> Vmin_M = -mv;
            Matrix<double> Vmax_M = mv;

            Vector<double> w = Vector<double>.Build.Dense(me, i => 0.9 - (0.7 * i / me));
            double c1 = 0.8;
            double c2 = 1.49;

            Random random = new Random();
            Func<double> randomFunc = () => random.NextDouble();

            Matrix<double> pos = Rmin_M + (Rmax_M - Rmin_M).PointwiseMultiply(DenseMatrix.Create(ps, D, (i, j) => randomFunc()));
            Matrix<double> vel = Vmin_M + (Vmax_M - Vmin_M).PointwiseMultiply(DenseMatrix.Create(ps, D, (i, j) => randomFunc()));

            List<Particle> particles = new List<Particle>();
            for (int i = 0; i < ps; ++i)
            {
                Particle p = new Particle
                {
                    Position = pos.Row(i),
                    Velocity = vel.Row(i),
                    PBest = pos.Row(i).Clone()
                };

                Vector<double> tri_var = pos.Row(i).SubVector(0, 6);
                double tA = tri_var[0];
                double tB = tri_var[1];
                double tC = tri_var[2];
                double deltav1 = tri_var[3];
                double alpha = tri_var[4];
                double beta = tri_var[5];

                TriImpRendez(tA, tB, tC, deltav1, alpha, beta, RT0, VT0, RP0, VP0, t0, out double cnt_fitness);

                p.Fitness = cnt_fitness;
                p.PBestFitness = cnt_fitness;

                particles.Add(p);
            }

            double gbestval = particles[0].PBestFitness;
            int minIndex = 0;
            for (int i = 1; i < ps; ++i)
            {
                if (particles[i].PBestFitness < gbestval)
                {
                    gbestval = particles[i].PBestFitness;
                    minIndex = i;
                }
            }

            Vector<double> gbest = particles[minIndex].PBest.Clone();
            Vector<double> real_best = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            for (int i = 1; i < me; ++i)
            {
                for (int k = 0; k < ps; ++k)
                {
                    particles[k].Velocity = w[i] * particles[k].Velocity +
                        c1 * (particles[k].PBest - particles[k].Position).PointwiseMultiply(DenseVector.Create(D, j => randomFunc())) +
                        c2 * (gbest - particles[k].Position).PointwiseMultiply(DenseVector.Create(D, j => randomFunc()));

                    particles[k].Velocity = particles[k].Velocity.PointwiseMaximum(-mv.Row(k)).PointwiseMinimum(mv.Row(k));
                    particles[k].Position = particles[k].Position + particles[k].Velocity;
                    particles[k].Position = particles[k].Position.PointwiseMaximum(Rmin_M.Row(k)).PointwiseMinimum(Rmax_M.Row(k));

                    Vector<double> tri_var2 = particles[k].Position.SubVector(0, 6);
                    double tA = tri_var2[0];
                    double tB = tri_var2[1];
                    double tC = tri_var2[2];
                    double deltav1 = tri_var2[3];
                    double alpha = tri_var2[4];
                    double beta = tri_var2[5];

                    TriImpRendez(tA, tB, tC, deltav1, alpha, beta, RT0, VT0, RP0, VP0, t0, out double cnt_fitness);
                    //Console.WriteLine($"dv={deltav1}");
                    //Console.WriteLine($"dv={cnt_fitness}");
                    particles[k].Fitness = cnt_fitness;
                    real_best = tri_var2;

                    if (particles[k].Fitness <= particles[k].PBestFitness)
                    {   
                        particles[k].PBest = particles[k].Position.Clone();
                        particles[k].PBestFitness = particles[k].Fitness;
                    }

                    if (particles[k].PBestFitness < gbestval)
                    {
                        gbest = particles[k].PBest.Clone();
                        gbestval = particles[k].PBestFitness;
                    }
                }
            }


            Vector<double> tri_vb = gbest;
            double tA1 = tri_vb[0];
            double tB1 = tri_vb[1];
            double tC1 = tri_vb[2];
            double deltav11 = tri_vb[3];
            double alpha1 = tri_vb[4];
            double beta1 = tri_vb[5];

            SolveImp(tA1, tB1, tC1, deltav11, alpha1, beta1, RT0, VT0, RP0, VP0, t0,
                out Vector<double> DeltaV1, out Vector<double> DeltaV2, out Vector<double> DeltaV3);
            //gbest.ToRowMatrix()


            //----------------------------------以下为原文件里的输出------------------------------------------------------------
            List<double> App_epoch_dv = new List<double> { tA1, tB1, tC1 };
            App_epoch_dv.Sort();
            double T_app_opt = App_epoch_dv[2];
            double Dv_sum = gbestval;

            int length = DeltaV1.Count; // 假设每个向量的长度相同
            // 初始化List存储分量
            List<double> App_Mag_dv_x = new List<double>();
            List<double> App_Mag_dv_y = new List<double>();
            List<double> App_Mag_dv_z = new List<double>();

            // 将分量填入对应的数组
            App_Mag_dv_x.Add(DeltaV1[0]);
            App_Mag_dv_x.Add(DeltaV2[0]);
            App_Mag_dv_x.Add(DeltaV3[0]);

            App_Mag_dv_y.Add(DeltaV1[1]);
            App_Mag_dv_y.Add(DeltaV2[1]);
            App_Mag_dv_y.Add(DeltaV3[1]);

            App_Mag_dv_z.Add(DeltaV1[2]);
            App_Mag_dv_z.Add(DeltaV2[2]);
            App_Mag_dv_z.Add(DeltaV3[2]);

            // 输出结果
            Console.WriteLine("App_epoch_dv: " + string.Join(", ", App_epoch_dv));
            Console.WriteLine("T_app_opt: " + string.Join(", ", T_app_opt));
            Console.WriteLine("Dv_sum: " + string.Join(", ", Dv_sum));
            Console.WriteLine("App_Mag_dv_x: " + string.Join(", ", App_Mag_dv_x));
            Console.WriteLine("App_Mag_dv_y: " + string.Join(", ", App_Mag_dv_y));
            Console.WriteLine("App_Mag_dv_z: " + string.Join(", ", App_Mag_dv_z));
            Console.ReadKey();

        }
    }
}