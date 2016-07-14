using System;
using System.Collections.Generic;
using System.Linq;

namespace NonclassicalVolterraIntegralEquations
{
    public class NonclassicalVolterraIntegralEquations
    {
        /// <summary>
        /// Implimentation of delegate lambda function for integral equation kernel: K(x,s).
        /// </summary>
        /// <param name="x">First variable.</param>
        /// <param name="s">Second variable.</param>
        /// <returns>Returns double value of kernel in x, s points.</returns>
        public delegate double Kernel(double x, double s);

        /// <summary>
        /// Implementation of delegate lambda function for known functions of integral equation: a(x), f(x), etc.
        /// </summary>
        /// <param name="x">Variable.</param>
        /// <returns>Returns double value of function in x point.</returns>
        public delegate double Function(double x);

        /// <summary>
        /// Implementation of right rectangle method algorithm for solving nonclassical Volterra integral equationof the first kind.
        /// </summary>
        /// <param name="kernel">Kernel of integral equation.</param>
        /// <param name="f">Right side function.</param>
        /// <param name="a">Lower limit function.</param>
        /// <param name="tmin">Left border of numeric grid.</param>
        /// <param name="tmax">Right border of numeric grid.</param>
        /// <param name="h">Step of numeric grid.</param>
        /// <param name="grid">Output array of grid's points.</param>
        /// <param name="solution">Output array of integral equation solution.</param>
        public static void NVEFirstKindRightRectangleMethod(Kernel kernel, Function f, Function a, double tmin, double tmax, double h, out double[] grid, out double[] solution)
        {
            int n = (int)Math.Round((tmax - tmin) / h, 0);

            double[] t = new double[n + 1];
            double[,] kernmatrix = new double[n + 1, n + 1];
            double[] rightside = new double[n + 1];

            grid = new double[n + 1];
            solution = new double[n + 1];

            int numRows = kernmatrix.GetLength(0), numCols = kernmatrix.GetLength(1);

            double[,] matrixforsolver = new double[numRows, numCols + 2];

            const double tiny = 0.00001;
            t[0] = tmin;
            t[n] = tmax;
            for (int i = 1; i < n; i++)
            {
                t[i] = t[i - 1] + h;
            }

            for (int i = 0; i <= n; i++)
            {
                kernmatrix[i, (int)Math.Round(a(t[i]) / h, 0) + 1] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t[(int)Math.Round(a(t[i]) / h, 0) + 1]);
                for (int k = (int)a(t[i] / h) + 2; k < i; k++)
                {
                    kernmatrix[i, k] += kernel(t[i], t[k]) * h;
                }
                kernmatrix[i, i] += kernel(t[i], t[i]) * h;
                rightside[i] += f(t[i]);
            }

            for (int i = 0; i <= n; i++)
            {
                for (int k = 0; k <= n; k++)
                {
                    matrixforsolver[i, k] = kernmatrix[i, k];
                }
                matrixforsolver[i, n + 1] = rightside[i];
            }

            for (int r = 0; r < numRows - 1; r++)
            {
                if (Math.Abs(matrixforsolver[r, r]) < tiny)
                {
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        if (Math.Abs(matrixforsolver[r2, r]) > tiny)
                        {
                            for (int c = 0; c <= numCols; c++)
                            {
                                double tmp = matrixforsolver[r, c];
                                matrixforsolver[r, c] = matrixforsolver[r2, c];
                                matrixforsolver[r2, c] = tmp;
                            }
                            break;
                        }
                    }
                }

                if (Math.Abs(matrixforsolver[r, r]) > tiny)
                {
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        double factor = -matrixforsolver[r2, r] / matrixforsolver[r, r];
                        for (int c = r; c <= numCols; c++)
                        {
                            matrixforsolver[r2, c] = matrixforsolver[r2, c] + factor * matrixforsolver[r, c];
                        }
                    }
                }
            }

            if (matrixforsolver[numRows - 1, numCols - 1] == 0)
            {
                bool allZeros = true;
                for (int c = 0; c <= numCols + 1; c++)
                {
                    if (matrixforsolver[numRows - 1, c] != 0)
                    {
                        allZeros = false;
                        break;
                    }
                }
                if (allZeros)
                {
                    throw new Exception("The solution is not unique.");
                }
                else
                {
                    throw new Exception("There is no solution.");
                }
            }
            else
            {
                for (int r = numRows - 1; r >= 0; r--)
                {
                    double tmp = matrixforsolver[r, numCols];
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        tmp -= matrixforsolver[r, r2] * matrixforsolver[r2, numCols + 1];
                    }
                    matrixforsolver[r, numCols + 1] = tmp / matrixforsolver[r, r];
                }

                for (int r = 0; r < numRows; r++)
                {
                    grid[r] = t[r];
                    solution[r] = matrixforsolver[r, numCols + 1];
                }
            }
        }

        /// <summary>
        /// Implementation of middle rectangle method algorithm for solving nonclassical Volterra integral equation of the first kind.
        /// </summary>
        /// <param name="kernel">Kernel of integral equation.</param>
        /// <param name="f">Right side function.</param>
        /// <param name="a">Lower limit function.</param>
        /// <param name="tmin">Left border of numeric grid.</param>
        /// <param name="tmax">Right border of numeric grid.</param>
        /// <param name="h">Step of numeric grid.</param>
        /// <param name="grid">Output array of grid's points.</param>
        /// <param name="solution">Output array of integral equation solution.</param>
        public static void NVEFirstKindMiddleRectangleMethod(Kernel kernel, Function f, Function a, double tmin, double tmax, double h, out double[] grid, out double[] solution)
        {
            int n = (int)Math.Round((tmax - tmin) / h, 0);

            double[] t = new double[n + 1];
            double[] t12 = new double[n + 1];
            double[] t12Sup = new double[n + 1];
            double[,] kernmatrix = new double[n + 1, n + 1];
            double[] rightside = new double[n + 1];

            double[] gridtemp = new double[n + 1];
            double[] solutiontemp = new double[n + 1];

            int numRows = kernmatrix.GetLength(0), numCols = kernmatrix.GetLength(1);

            double[,] matrixforsolver = new double[numRows, numCols + 2];

            const double tiny = 0.00001;

            t[0] = t12[0] = t12Sup[0] = tmin;
            t[n] = tmax;

            for (int i = 1; i < n; i++)
            {
                t[i] = t[i - 1] + h;
            }
            for (int i = 1; i <= n; i++)
            {
                t12[i] = t[0] + (i - 0.5) * h;
                t12Sup[i] = (t[(int)Math.Round(a(t[i]) / h) + 1] + a(t[i])) / 2;
            }
            List<int> toCalculate = new List<int>();
            List<Tuple<double, double>> output = new List<Tuple<double, double>>();

            for (int i = 0; i <= n; i++)
            {
                if (t12.Contains(t12Sup[i]))
                    kernmatrix[i, Array.IndexOf(t12, t12Sup[i])] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]);
                else
                {
                    if (t12Sup[i] < t12[0])
                    {
                        kernmatrix[i, 0] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                            (1 + (t12Sup[i] - t[0]) / (t[1] - t[0]));
                        kernmatrix[i, 1] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                            (-(t12Sup[i] - t[0]) / (t[1] - t[0]));
                    }
                    if (t12Sup[i] > t12[n])
                    {
                        kernmatrix[i, n - 1] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                            (1 + (t12Sup[i] - t[n - 1]) / (t[n] - t[n - 1]));
                        kernmatrix[i, n] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                            (-(t12Sup[i] - t[n - 1]) / (t[n] - t[n - 1]));
                    }
                    for (int j = 0; j < n; j++)
                    {
                        if (t12Sup[j] > t12[j] && t12Sup[j] < t12[j + 1])
                        {
                            kernmatrix[i, j] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                            (1 + (t12Sup[i] - t[j]) / (t[j] - t[j + 1]));
                            kernmatrix[i, j + 1] += (t[(int)Math.Round(a(t[i]) / h, 0) + 1] - a(t[i])) * kernel(t[i], t12Sup[i]) *
                                                (-(t12Sup[i] - t[j]) / (t[j] - t[j + 1]));
                        }
                    }
                    toCalculate.Add(i);
                }
                for (int k = (int)a(t[i] / h) + 2; k < i; k++)
                {
                    kernmatrix[i, k] += kernel(t[i], t12[k]) * h;
                }
                kernmatrix[i, i] += kernel(t[i], t12[i]) * h;
                rightside[i] += f(t[i]);

            }
            for (int i = 0; i <= n; i++)
            {
                for (int k = 0; k <= n; k++)
                {
                    matrixforsolver[i, k] = kernmatrix[i, k];
                }
                matrixforsolver[i, n + 1] = rightside[i];
            }

            for (int r = 0; r < numRows - 1; r++)
            {
                if (Math.Abs(matrixforsolver[r, r]) < tiny)
                {
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        if (Math.Abs(matrixforsolver[r2, r]) > tiny)
                        {
                            for (int c = 0; c <= numCols; c++)
                            {
                                double tmp = matrixforsolver[r, c];
                                matrixforsolver[r, c] = matrixforsolver[r2, c];
                                matrixforsolver[r2, c] = tmp;
                            }
                            break;
                        }
                    }
                }

                if (Math.Abs(matrixforsolver[r, r]) > tiny)
                {
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        double factor = -matrixforsolver[r2, r] / matrixforsolver[r, r];
                        for (int c = r; c <= numCols; c++)
                        {
                            matrixforsolver[r2, c] = matrixforsolver[r2, c] + factor * matrixforsolver[r, c];
                        }
                    }
                }
            }

            if (matrixforsolver[numRows - 1, numCols - 1] == 0)
            {
                bool allZeros = true;
                for (int c = 0; c <= numCols + 1; c++)
                {
                    if (matrixforsolver[numRows - 1, c] != 0)
                    {
                        allZeros = false;
                        break;
                    }
                }
                if (allZeros)
                {
                    throw new Exception("The solution is not unique.");
                }
                else
                {
                    throw new Exception("There is no solution.");
                }
            }
            else
            {
                for (int r = numRows - 1; r >= 0; r--)
                {
                    double tmp = matrixforsolver[r, numCols];
                    for (int r2 = r + 1; r2 < numRows; r2++)
                    {
                        tmp -= matrixforsolver[r, r2] * matrixforsolver[r2, numCols + 1];
                    }
                    matrixforsolver[r, numCols + 1] = tmp / matrixforsolver[r, r];
                }

                for (int r = 0; r < numRows; r++)
                {
                    output.Add(new Tuple<double, double>(t[r], matrixforsolver[r, numCols + 1]));
                    gridtemp[r] = t[r];
                    solutiontemp[r] = matrixforsolver[r, numCols + 1];
                }
            }

            for (int i = 0; i <= n; i++)
            {
                output.Add(new Tuple<double, double>(gridtemp[i], solutiontemp[i]));
            }
            
            foreach (var i in toCalculate)
            {
                for (int j = 0; j < toCalculate.Count; j++)
                {
                    if (t12Sup[i] > t12[j] && t12Sup[i] < t12[j + 1])
                    {
                        output.Add(new Tuple<double, double>(t12Sup[i], solutiontemp[j] + (solutiontemp[j + 1] - solutiontemp[j]) * (t12Sup[i] - t12[j]) / (t12[j + 1] - t12[j])));
                    }
                }
            }

            output.Sort((y, x) => y.Item1.CompareTo(x.Item1));
            output = output.Distinct().ToList();

            grid = new double[output.Count];
            solution = new double[output.Count];

            for (int i = 0; i < output.Count; i++)
            {
                grid[i] = output[i].Item1;
                solution[i] = output[i].Item2;
            }
        }
    }
}
