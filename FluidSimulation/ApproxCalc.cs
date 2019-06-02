using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FluidSimulation
{
    class ApproxCalc
    {
        







        public double[,] CalcPressure(double[,] data)
        {
            return null;
        }

        public double[,] Gradient2D(double[,] input, double[,] points, double nzero, double radius)
        {
            var output = new double[input.GetLength(0), 2];
            for (int i = 0; i < input.GetLength(0); i++)
            {
                var weights = Weights(points, i, radius);
                for (int j = 0; j < input.GetLength(0); j++)
                {
                    var temp = weights[j] * ()
                }
            }
        }

        public double[] Weights(double[,] input, int basepointindex, double radius)
        {
            var output = new double[input.GetLength(0)];
            var x = input[basepointindex, 0];
            var y = input[basepointindex, 1];
            for (int i = 0; i < basepointindex; i++)
            {
                var dis = Distance(x, y, input[i, 0], input[i, 1]);
                output[i] = Math.Max(0, radius / dis - 1);
            }
            for (int i = basepointindex + 1; i < input.GetLength(0); i++)
            {
                var dis = Distance(x, y, input[i, 0], input[i, 1]);
                output[i] = Math.Max(0, radius / dis - 1);
            }
            return output;
        }

        public double[] 

        private double Nzero(double[,] input, int basepointindex, double radius)
        {
            var weights = Weights(input, basepointindex, radius);
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
    }
}
