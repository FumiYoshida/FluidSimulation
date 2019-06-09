using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Alea;
using Alea.Parallel;

namespace FluidSimulation
{
    class ApproxCalc
    {
        // ここから自動生成される値
        private double nzero;
        private double lambdazero;
        private double radius;
        private double timespan;
        private double gravity;
        private double vu;
        private double rho;
        public double[][,] data;
        private int[] gridtoindexes;
        private int[,,] gridsupporter;
        public bool[] iswall;
        public bool[] isdummy;
        private Gpu gpu;
        // ここまで自動生成される値

        // デバッグ用変数
        public double[] pressure;

        /*
        public ApproxCalc(double[,] points, int pointindex, double radius)
        {
            this.radius = radius;
            nzero = Nzero(points, pointindex);
            lambdazero = Lambdazero(points, pointindex);
            timespan = 0.0001; // 1タイムステップの時間(s)
            gravity = 0; // 重力加速度(m / s^2)
            vu = 1;// 動粘度(m^2 / s)
            rho = 1; // 密度(kg / m^3)
        }
        */

        public ApproxCalc()
        {
            Pillar();
            radius = 0.31;
            int pointindex = 280;
            nzero = Nzero(data[0], pointindex);
            lambdazero = Lambdazero(data[0], pointindex);
            timespan = 0.003; // 1タイムステップの時間(s)
            gravity = 9.8; // 重力加速度(m / s^2)
            vu = 1;// 動粘度(m^2 / s)
            rho = 1; // 密度(kg / m^3)
            gpu = Gpu.Default;
        }

        private void PutintoGrid()
        {
            int pointnum = data[0].GetLength(0);
            var xs = new double[pointnum];
            var ys = new double[pointnum];
            for (int i = 0; i < pointnum; i++)
            {
                xs[i] = data[0][i, 0];
                ys[i] = data[0][i, 0];
            }
            Array.Sort(xs);
            Array.Sort(ys);
            int PlacetoIndex(double x) => (int)Math.Max(0, Math.Floor(2 + (x / radius)));
            int yoko = PlacetoIndex(xs[pointnum - 1]) + 1;
            int tate = PlacetoIndex(ys[pointnum - 1]) + 1;
            var tempgrid = new List<int>[yoko, tate];
            for (int i = 0; i < yoko; i++)
            {
                for (int j = 0; j < tate; j++)
                {
                    tempgrid[i, j] = new List<int>();
                }
            }
            for (int i = 0; i < pointnum; i++)
            {
                int xindex = (int)Math.Max(0, Math.Floor(2 + (data[0][i, 0] / radius)));
                int yindex = (int)Math.Max(0, Math.Floor(2 + (data[0][i, 1] / radius)));
                for (int j = xindex - 1; j < xindex + 2; j++)
                {
                    if (j > -1 && j < yoko)
                    {
                        for (int k = yindex - 1; k < yindex + 2; k++)
                        {
                            if (k > -1 && k < tate)
                            {
                                tempgrid[j, k].Add(i);
                            }
                        }
                    }
                }
            }
            // 各グリッドに対応するgridtoindexのindexの開始地点がgridsupporter[i, j, 0]に、
            // そのグリッド内の粒子の数がgridsupporter[i, j, 1] に入る。
            gridsupporter = new int[yoko, tate, 2];
            var tempgrid2 = new List<int>();
            int serial = 0;
            for (int i = 0; i < yoko; i++)
            {
                for (int j = 0; j < tate; j++)
                {
                    gridsupporter[i, j, 0] = serial;
                    gridsupporter[i, j, 1] = tempgrid[i, j].Count();
                    serial += gridsupporter[i, j, 1];
                    tempgrid2.AddRange(tempgrid[i, j]);
                }
            }
            gridtoindexes = tempgrid2.ToArray();
        }

        public double[,] Gradient2D(double[] input, double[,] points)
        {
            // inputのgradientを取る。
            // points[i, 0], points[i, 1] はそれぞれinput[i, j] に対応するx座標、y座標
            // output[i, 0], output[i, 1] はそれぞれinputのgradientのx方向成分、y方向成分
            int pointnum = points.GetLength(0);
            var output = new double[pointnum, 2];
            var weights = new double[pointnum, pointnum];
            var radius = this.radius; // Aleaではクラスの変数を直接扱えない
            var nzero = this.nzero;
            var gridtoindexes = this.gridtoindexes;
            var gridsupporter = this.gridsupporter;
            gpu.For(0, input.Length, i =>
            {
                var x = points[i, 0];
                var y = points[i, 1];
                var gridx = (int)Math.Max(0, Math.Floor(2 + (x / radius)));
                var gridy = (int)Math.Max(0, Math.Floor(2 + (y / radius)));
                double tempx = 0;
                double tempy = 0;
                int fromindex = gridsupporter[gridx, gridy, 0];
                int pnum = gridsupporter[gridx, gridy, 1];
                for (int j0 = 0; j0 < pnum; j0++)
                {
                    int j = gridtoindexes[fromindex + j0];
                    if (i != j)
                    {
                        var xsa = x - points[j, 0];
                        var ysa = y - points[j, 1];
                        var dis = Math.Sqrt(xsa * xsa + ysa * ysa);
                        var weight = Math.Max(0, radius / dis - 1);
                        var xvec = points[j, 0] - points[i, 0];
                        var yvec = points[j, 1] - points[i, 1];
                        var temp = weight * (input[j] - input[i]) / (xvec * xvec + yvec * yvec);
                        tempx += temp * xvec;
                        tempy += temp * yvec;
                    }
                }
                output[i, 0] = tempx * 2 / nzero;
                output[i, 1] = tempy * 2 / nzero;
            });
            return output;
        }

        public double[,] Gradient2DforPressure(double[] input, double[,] points)
        {
            // inputのgradientを取る。
            // points[i, 0], points[i, 1] はそれぞれinput[i, j] に対応するx座標、y座標
            // output[i, 0], output[i, 1] はそれぞれinputのgradientのx方向成分、y方向成分
            int pointnum = points.GetLength(0);
            var output = new double[pointnum, 2];
            var weights = new double[pointnum, pointnum];

            var radius = this.radius; // Aleaではクラスの変数を直接扱えない
            var nzero = this.nzero;
            var gridtoindexes = this.gridtoindexes;
            var gridsupporter = this.gridsupporter;
            gpu.For(0, input.Length, i =>
            {
                var x = points[i, 0];
                var y = points[i, 1];
                var gridx = (int)Math.Max(0, Math.Floor(2 + (x / radius)));
                var gridy = (int)Math.Max(0, Math.Floor(2 + (y / radius)));
                double tempx = 0;
                double tempy = 0;
                double mininput = input[i];
                int fromindex = gridsupporter[gridx, gridy, 0];
                int pnum = gridsupporter[gridx, gridy, 1];
                for (int j0 = 0; j0 < pnum; j0++)
                {
                    int j = gridtoindexes[fromindex + j0];
                    mininput = Math.Min(input[j], mininput);
                }
                for (int j0 = 0; j0 < pnum; j0++)
                {
                    int j = gridtoindexes[fromindex + j0];
                    if (i != j)
                    {
                        var xsa = x - points[j, 0];
                        var ysa = y - points[j, 1];
                        var dis = Math.Sqrt(xsa * xsa + ysa * ysa);
                        var weight = Math.Max(0, radius / dis - 1);
                        var xvec = points[j, 0] - points[i, 0];
                        var yvec = points[j, 1] - points[i, 1];
                        var temp = weight * (input[j] - mininput) / (xvec * xvec + yvec * yvec);
                        tempx += temp * xvec;
                        tempy += temp * yvec;
                    }
                }
                output[i, 0] = tempx * 2 / nzero;
                output[i, 1] = tempy * 2 / nzero;
            });
            return output;
        }

        private double[] Weights(double[,] points, int basepointindex)
        {
            var output = new double[points.GetLength(0)];
            var x = points[basepointindex, 0];
            var y = points[basepointindex, 1];
            for (int i = 0; i < basepointindex; i++)
            {
                var dis = Distance(x, y, points[i, 0], points[i, 1]);
                output[i] = Math.Min(10, Math.Max(0, radius / dis - 1));
            }
            for (int i = basepointindex + 1; i < points.GetLength(0); i++)
            {
                var dis = Distance(x, y, points[i, 0], points[i, 1]);
                output[i] = Math.Max(0, radius / dis - 1);
            }
            return output;
        }

        private double Nzero(double[,] points, int basepointindex)
        {
            var weights = Weights(points, basepointindex);
            double sumweights = 0;
            for (int i = 0; i < weights.GetLength(0); i++)
            {
                sumweights += weights[i];
            }
            return sumweights;
        }

        private double Distance(double x1, double y1, double x2, double y2)
        {
            var xdif = x1 - x2;
            var ydif = y1 - y2;
            return Math.Sqrt(xdif * xdif + ydif * ydif);
        }

        public double Lambdazero(double[,] points, int basepointindex)
        {
            int pointnum = points.GetLength(0);
            var weights = Weights(points, basepointindex);
            double sumweights = 0;
            for (int i = 0; i < pointnum; i++)
            {
                sumweights += weights[i];
            }
            double output = 0;
            for (int i = 0; i < pointnum; i++)
            {
                double dis = Distance(points[basepointindex, 0], points[basepointindex, 1], points[i, 0], points[i, 1]);
                output += dis * dis * weights[i];
            }
            output /= sumweights;
            return output;
        }

        public double[,] Laplace2D(double[,] input, double[,] points)
        {
            int pointnum = points.GetLength(0);
            var output = new double[pointnum, 2];
            var nzero = this.nzero; // Aleaではクラス変数を直接扱えない
            var lambdazero = this.lambdazero;
            var radius = this.radius;
            var weights = new double[pointnum, pointnum];
            var gridtoindexes = this.gridtoindexes;
            var gridsupporter = this.gridsupporter;
            var nzlmz = 4 / nzero / lambdazero;
            gpu.For(0, pointnum, i => {
                var x = points[i, 0];
                var y = points[i, 1];
                var gridx = (int)Math.Max(0, Math.Floor(2 + (x / radius)));
                var gridy = (int)Math.Max(0, Math.Floor(2 + (y / radius)));
                double tempx = 0;
                double tempy = 0;
                int fromindex = gridsupporter[gridx, gridy, 0];
                int pnum = gridsupporter[gridx, gridy, 1];
                for (int j0 = 0; j0 < pnum; j0++)
                {
                    int j = gridtoindexes[fromindex + j0];
                    if (i != j)
                    {
                        var xsa = x - points[j, 0];
                        var ysa = y - points[j, 1];
                        var dis = Math.Sqrt(xsa * xsa + ysa * ysa);
                        var weight = Math.Max(0, radius / dis - 1);
                        tempx += (input[j, 0] - input[i, 0]) * weight;
                        tempy += (input[j, 1] - input[i, 1]) * weight;
                    }
                }
                output[i, 0] = tempx * 4 / nzero / lambdazero;
                output[i, 1] = tempy * 4 / nzero / lambdazero;
            });
            return output;
        }

        public void Pillar()
        {
            int pointnum = 1672;
            var points = new double[pointnum, 2];
            iswall = new bool[pointnum];
            isdummy = new bool[pointnum];
            int serial = 0;
            for (int i = 0; i < 10; i++)
            {
                for (int j = 0; j < 50; j++)
                {
                    points[serial, 0] = i * 0.01;
                    points[serial, 1] = j * 0.01;
                    serial++;
                }
            }
            for (int i = 10; i < 60; i++)
            {
                for (int j = 0; j < 10; j++)
                {
                    points[serial, 0] = i * 0.01;
                    points[serial, 1] = j * 0.01;
                    serial++;
                }
            }

            // 両端に壁粒子とダミー粒子を配置
            for (int j = 0; j < 50; j++)
            {
                points[serial, 0] = -0.01;
                points[serial, 1] = j * 0.01;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = -0.02;
                points[serial, 1] = j * 0.01;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = 0.6;
                points[serial, 1] = j * 0.01;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = 0.61;
                points[serial, 1] = j * 0.01;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = -0.03;
                points[serial, 1] = j * 0.01;
                isdummy[serial] = true;
                serial++;
                points[serial, 0] = -0.04;
                points[serial, 1] = j * 0.01;
                isdummy[serial] = true;
                serial++;
                points[serial, 0] = 0.62;
                points[serial, 1] = j * 0.01;
                isdummy[serial] = true;
                serial++;
                points[serial, 0] = 0.63;
                points[serial, 1] = j * 0.01;
                isdummy[serial] = true;
                serial++;
            }

            // 下に壁粒子とダミー粒子を配置

            for (int i = -4; i < 64; i++)
            {
                points[serial, 0] = i * 0.01;
                points[serial, 1] = -0.01;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = i * 0.01;
                points[serial, 1] = -0.02;
                iswall[serial] = true;
                serial++;
                points[serial, 0] = i * 0.01;
                points[serial, 1] = -0.03;
                isdummy[serial] = true;
                serial++;
                points[serial, 0] = i * 0.01;
                points[serial, 1] = -0.04;
                isdummy[serial] = true;
                serial++;
            }

            data = new double[2][,] { points, new double[pointnum, 2] };
        }

        private double[,] Velocity2D(double[,] points, double[] pressure, double[,] velocity)
        {
            int pointnum = pressure.GetLength(0);
            var dudt = new double[pointnum, 2];
            var gradpressure = Gradient2D(pressure, points);
            for (int i = 0; i < pointnum; i++)
            {
                if (!(iswall[i] | isdummy[i]))
                {
                    dudt[i, 0] = -gradpressure[i, 0] / rho;
                    dudt[i, 1] = -gradpressure[i, 1] / rho;
                }
            }
            var output = new double[pointnum, 2];
            for (int i = 0; i < pointnum; i++)
            {
                output[i, 0] = velocity[i, 0] + timespan * dudt[i, 0];
                output[i, 1] = velocity[i, 1] + timespan * dudt[i, 1];
            }
            return output;
        }

        private double[,] Velocity2D(double[,] points, double[,] velocity)
        {
            int pointnum = points.GetLength(0);
            var dudt = new double[pointnum, 2];
            var laplacevelocity = Laplace2D(velocity, points);
            for (int i = 0; i < pointnum; i++)
            {
                if (!(iswall[i] | isdummy[i]))
                {
                    dudt[i, 0] += vu * laplacevelocity[i, 0];
                    dudt[i, 1] += vu * laplacevelocity[i, 1] - gravity;
                }
            }
            var output = new double[pointnum, 2];
            for (int i = 0; i < pointnum; i++)
            {
                output[i, 0] = velocity[i, 0] + timespan * dudt[i, 0];
                output[i, 1] = velocity[i, 1] + timespan * dudt[i, 1];
            }
            return output;
        }

        private double[,] MovePoints2D(double[,] points, double[,] velocity)
        {
            var output = new double[points.GetLength(0), 2];
            for (int i = 0; i < points.GetLength(0); i++)
            {
                output[i, 0] = points[i, 0] + timespan * velocity[i, 0];
                output[i, 1] = points[i, 1] + timespan * velocity[i, 1];
            }
            return output;
        }

        public void Next()
        {
            var points = data[0];
            var velocity = data[1];

            PutintoGrid();
            // まず圧力を考えずに粒子を移動させる
            var tempvelocity = Velocity2D(points, velocity);
            var temppoints = MovePoints2D(points, tempvelocity);

            PutintoGrid();
            // 移動後の粒子の位置をもとに圧力を計算する
            var pressure = Pressure(temppoints);
            this.pressure = pressure;

            // 圧力の分粒子を速度を更新し、最終的な速度を計算する
            var nextvelocity = Velocity2D(temppoints, pressure, tempvelocity);

            // 速度をもとに粒子の位置を移動させる
            var nextpoints = MovePoints2D(points, nextvelocity);

            data = new double[2][,] { nextpoints, nextvelocity };
        }

        private double[] Pressure(double[,] points)
        {
            // まず粒子数密度を求める
            var pointnum = points.GetLength(0);
            var parden = new double[pointnum];
            var issurface = new bool[pointnum];
            double gamma = 0.2; // これをかけて安定させるらしい
            var nzero = this.nzero; // Aleaではクラス変数を直接扱えない
            var lambdazero = this.lambdazero;
            var radius = this.radius;
            var weights = new double[pointnum, pointnum];
            var multi = gamma / nzero / timespan / timespan;
            var isnotdummyorwall = new bool[pointnum];
            for (int i = 0; i < pointnum; i++)
            {
                isnotdummyorwall[i] = !(isdummy[i] | iswall[i]);
            }
            gpu.For(0, pointnum, i => {
                var x = points[i, 0];
                var y = points[i, 1];
                double tempparden = 0;
                for (int j = 0; j < pointnum; j++)
                {
                    if (i != j)
                    {
                        var xsa = x - points[j, 0];
                        var ysa = y - points[j, 1];
                        var dis = Math.Sqrt(xsa * xsa + ysa * ysa);
                        weights[i, j] = Math.Max(0, radius / dis - 1);
                        tempparden += weights[i, j];
                    }
                }
                issurface[i] = isnotdummyorwall[i] & tempparden < nzero * 0.97;
                parden[i] = multi * (tempparden - nzero);
            });
            
            for (int j = 0; j < pointnum; j++)
            {
                if (isdummy[j] | issurface[j])
                {
                    for (int i = 0; i < pointnum; i++)
                    {
                        weights[i, j] = 0; // 圧力を持たないダミー粒子と、表面(=圧力0)の粒子は連立方程式に関わらない
                    }
                }
            }
            // weights * pressure = parden となる。weightsは疎行列なので共役勾配法(conjugate gradient method) を用いてpressureを求める。
            var pressure = ConjugateGradientMethod(weights, parden);
            for (int i = 0; i < pointnum; i++)
            {
                if (isdummy[i] | issurface[i])
                {
                    pressure[i] = 0;
                }
                else if (pressure[i] < 0)
                {
                    pressure[i] = 0;
                }
            }
            return pressure;
        }

        private double[] ConjugateGradientMethod(double[,] weight, double[] product)
        {
            int tate = product.Length;
            int tako = weight.GetLength(0);
            var output = new double[tako];
            var rest = new double[tako];
            Array.Copy(product, rest, tako);
            var p = new double[tako];
            Array.Copy(product, p, tako);
            double threshold = rest.Average() / 100;
            while (rest.Average() > threshold)
            {
                double[] wp = Product(weight, p);
                double a = Product(rest, p) / Product(p, wp);
                output = Add(output, Product(a, p));
                double[] nextrest = Subtract(rest, Product(a, wp));
                double b = Product(nextrest, nextrest) / Product(rest, rest);
                p = Add(nextrest, Product(b, p));
                Array.Copy(nextrest, rest, rest.Length);
            }
            return output;
        }


        // ここから一般的な関数
        private double Product(double[] x, double[] y)
        {
            double output = 0;
            for (int i = 0; i < Math.Min(x.Length, y.Length); i++)
            {
                output += x[i] * y[i];
            }
            return output;
        }

        private double[] Product(double[,] x, double[] y)
        {
            int tate = x.GetLength(0);
            int tako = Math.Min(y.Length, x.GetLength(1));
            var output = new double[tate];
            gpu.For(0, tate, i => {
                double temp = 0;
                for (int j = 0; j < tako; j++)
                {
                    temp += x[i, j] * y[j];
                }
                output[i] = temp;
            });
            return output;
        }

        private double[] Product(double x, double[] y)
        {
            var output = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                output[i] = x * y[i];
            }
            return output;
        }

        private double[] Subtract(double[] x, double[] y)
        {
            int tate = Math.Min(x.Length, y.Length);
            var output = new double[tate];
            for (int i = 0; i < tate; i++)
            {
                output[i] = x[i] - y[i];
            }
            return output;
        }

        private double[] Add(double[] x, double[] y)
        {
            int tate = Math.Min(x.Length, y.Length);
            var output = new double[tate];
            for (int i = 0; i < tate; i++)
            {
                output[i] = x[i] + y[i];
            }
            return output;
        }

        // ここまで一般的な関数
    }
}
