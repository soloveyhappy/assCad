#region Usings

using System;
using System.Collections.Generic;
using TimeZYX.Model.Spatial;
using System.Threading;
using System.Text;

#endregion

namespace TimeZYX.Model.Utils
{

    public static class MathUtils
    {
        public static System.Random RndObject = new System.Random(System.DateTime.Now.Millisecond);

        public const double MachineEpsilon = 5E-16;
        public const double MaxRealNumber = 1E300;
        public const double MinRealNumber = 1E-300;
        public static readonly float PI_180F = Convert.ToSingle(Math.PI / 180d);
        public static readonly double PI_180D = Math.PI / 180d;

        public static double RandomReal()
        {
            double r = 0;
            lock (RndObject) { r = RndObject.NextDouble(); }
            return r;
        }
        public static int RandomInteger(int N)
        {
            int r = 0;
            lock (RndObject) { r = RndObject.Next(N); }
            return r;
        }
        public static double Sqr(double X)
        {
            return X * X;
        }
        public static void CubicRoots(double a, double b, double c, double d, out double[] troot)
        {
            double des;
            List<double> troots = new List<double>(3);
            if (a == 0)
            {
                des = c * c - b * d * 4;
                if (des < 0)
                {
                    troot = new double[0];
                    return;
                }
                else
                {
                    troot = new double[] { (-c - Math.Sqrt(des)) / (2 * b), (-c + Math.Sqrt(des)) / (2 * b) };
                    return;
                }
            }
            double q = (3 * a * c - b * b) / (9 * a * a);
            double r = (9 * a * b * c - 27 * a * a * d - 2 * b * b * b) / (54 * a * a * a);
            des = q * q * q + r * r;
            double s, t;
            double x1, x2;
            if (des > 0)
            {
                s = Math.Sign(r + Math.Sqrt(des)) * Math.Pow(Math.Abs(r + Math.Sqrt(des)), 1.0 / 3);
                t = Math.Sign(r - Math.Sqrt(des)) * Math.Pow(Math.Abs(r - Math.Sqrt(des)), 1.0 / 3);
                troot = new double[] { s + t - b / 3 / a };
                return;
            }
            if (des == 0)
            {

                s = Math.Pow(-q * q * q, 1.0 / 6);
                x1 = 2 * s - b / 3 / a;
                x2 = -s - b / 3 / a;
                if (x1 == x2)
                    troot = new double[] { x1 };
                else troot = new double[] { x1, x2 };
                return;
            }

            s = Math.Sqrt(-q * q * q);
            t = Math.Acos(r / s) / 3;
            s = Math.Pow(s, 1.0 / 3);
            troot = new double[] { 2 * s * Math.Cos(t) - b / 3 / a, s * (-Math.Sqrt(3) * Math.Sin(t) - Math.Cos(t)) - b / 3 / a, s * (Math.Sqrt(3) * Math.Sin(t) - Math.Cos(t)) - b / 3 / a };
            return;
        }
        public static double AbsComplex(Complex z)
        {
            double w;
            double xabs;
            double yabs;
            double v;

            xabs = System.Math.Abs((double)((sbyte)z.x));
            yabs = System.Math.Abs((double)((sbyte)z.y));
            w = xabs > yabs ? xabs : yabs;
            v = xabs < yabs ? xabs : yabs;
            if (v == 0)
                return w;
            else
            {
                double t = v / w;
                return w * System.Math.Sqrt(1 + t * t);
            }
        }
        public static Complex Conj(Complex z)
        {
            return new Complex(z.x, -z.y);
        }
        public static Complex CSqr(Complex z)
        {
            return new Complex(z.x * z.x - z.y * z.y, 2 * z.x * z.y);
        }
        public static Complex CSqrt(Complex z)
        {
            double r = Math.Sqrt(z.x * z.x + z.y * z.y);
            if (r==0) return new Complex(0,0);
            double fi = Math.Acos(z.x / r);
            if (z.y < 0) fi = 2 * Math.PI - fi;
            return new Complex(Math.Sqrt(r) * Math.Cos(fi / 2), Math.Sqrt(r) * Math.Sin(fi / 2));
        }
        public static Complex CSqrtn(Complex z, int n)
        {
            double r = Math.Sqrt(z.x * z.x + z.y * z.y);
            if (r == 0) return new Complex(0, 0);
            double fi = Math.Acos(z.x / r);
            if (z.y < 0) fi = 2 * Math.PI - fi;
            return new Complex(Math.Pow(r, 1.0 / n) * Math.Cos(fi / n + 2 * Math.PI / n), Math.Pow(r, 1.0 / n) * Math.Sin(fi / n + 2 * Math.PI / n));
        }
        public static Complex[] CSqrt3(Complex z)
        {
            int n = 3;
            Complex[] x = new Complex[3];
            double r = Math.Sqrt(z.x * z.x + z.y * z.y);
            if (r == 0) return x;
            double fi = Math.Acos(z.x / r);
            if (z.y < 0) fi = 2 * Math.PI - fi;
            r = Math.Pow(r, 1.0 / n);
            fi /= n;
            for (int i = 0; i < n;i++ )
                x[i] = new Complex(r * Math.Cos(fi + i * 2 * Math.PI / n), r * Math.Sin(fi+i * 2 * Math.PI / n));
           
            return x;
        }
    }


    public static class MathFunctions
    {
        #region Возведение в квадрат

        public static int Pow2(int x)
        {
            return x*x;
        }

        public static float Pow2(float x)
        {
            return x*x;
        }

        public static double Pow2(double x)
        {
            return x*x;
        }

        #endregion

        #region Растояние между двумя точками...
        /// <summary>
        /// Растояние между двумя точками в 3D пространстве
        /// </summary>

        public static double Distance2Points(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            return Math.Sqrt(Pow2(x1 - x2) + Pow2(y1 - y2) + Pow2(z1 - z2));
        }
        /// <summary>
        /// Растояние между точками в 2d пространстве
        /// </summary>
        public static double Distance2Points(double x1, double y1, double x2, double y2)
        {
            return Math.Sqrt(Pow2(x1 - x2) + Pow2(y1 - y2));
        }
        #endregion 

        #region Матричные операции

        /// <summary>
        /// Конвертировать в квадратную матрицу
        /// </summary>
        /// <param name="matrix">Одномереный массив правильной размерности</param>
        /// <returns>Двумерная матрица</returns>
        public static float[,] ToMatrixForm(float[] matrix)
        {
            int n = (int) Math.Sqrt(matrix.Length);

            if (matrix.Length != Pow2(n))
                throw new ArgumentException("Non-square matrix argument");

            float[,] matrix2 = new float[n,n];

            int k = 0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    matrix2[i, j] = matrix[k++];

            return matrix2;
        }

        /// <summary>
        /// Транспонировать матрицу
        /// </summary>
        /// <param name="matrix">Матрица</param>
        /// <returns>Транспонированная матрица</returns>
        public static float[,] Transpose(float[,] matrix)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n)
                throw new ArgumentException("Non-square matrix argument");

            float[,] matrix2 = new float[n,n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    matrix2[j, i] = matrix[i, j];

            return matrix2;
        }

        public static float[] GaussSolve(float[] matrix, float[] right)
        {
            int n = right.Length;
            if (matrix.Length != Pow2(n))
                throw new ArgumentException("Incompatible arguments");

            return GaussSolve(ToMatrixForm(matrix), right);
        }

        public static float[] GaussSolve(float[,] matrix, float[] right)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n || right.Length != n)
                throw new ArgumentException("Incompatible arguments");

            float[] vector = (float[]) right.Clone();
            float[,] m = (float[,]) matrix.Clone();

            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++)
                    if (i != j)
                    {
                        float a = m[i, j]/m[j, j];
                        vector[i] -= vector[j]*a;
                        for (int k = j; k < n; k++)
                            m[i, k] -= m[j, k]*a;
                    }

            float[] result = new float[n];
            for (int i = 0; i < n; i++)
                result[i] = vector[i]/m[i, i];

            return result;
        }


        public static float[] MatrixMULVector(float[] matrix, float[] vector)
        {
            int n = vector.Length;
            if (matrix.Length != Pow2(n))
                throw new ArgumentException("Incompatible arguments");

            return MatrixMULVector(ToMatrixForm(matrix), vector);
        }

        public static float[] MatrixMULVector(float[,] matrix, float[] vector)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n || vector.Length != n)
                throw new ArgumentException("Incompatible arguments");

            float[] result = new float[n];

            for (int i = 0; i < n; i++)
            {
                result[i] = 0;
                for (int j = 0; j < n; j++)
                    result[i] += matrix[i, j]*vector[j];
            }

            return result;
        }

        /// <summary>
        /// Получить минимум/максимум/среднее по вектору значений
        /// </summary>
        /// <param name="vector">Вектор</param>
		/// <param name="min">Out-параметр: минимум. Если минимум найден не был, то float.MinValue </param>
        /// <param name="max">Out-параметр: максимум. Если максимум найден не был, то float.MaxValue</param>
        /// <param name="avg">Out-параметр: среднее значение. Если количество значений = 0, то возвращается NaN</param>
        public static void GetMinMaxAvg(IEnumerable<float> vector, out float min, out float max, out float avg)
        {
            min = float.MaxValue;
            max = float.MinValue;

            long counter = 1;
            Double sum = 0.0f;
			if (vector != null)
			{
				counter = 0;
				foreach (float v in vector)
				{
					if (v < min)
						min = v;
					if (v > max)
						max = v;
					sum += v;
					counter++;
				}
			}
            avg = (Single)(sum / counter);
        }

        #endregion

        #region Геометрия полигонов

        public static bool Between(float x, float a, float b)
        {
            return (a <= x && x <= b) || (b <= x && x <= a);
        }

        public static bool Between<T>(T x, T a, T b) where T : IComparable<T>
        {
            return (a.CompareTo(x) <= 0 && x.CompareTo(b) <= 0) || (b.CompareTo(x) <= 0 && x.CompareTo(a) <= 0);
        }

        public static bool PointInRegion2D(IVertexDouble v, IList<IVertexDouble> region)
        {
            if (region.Count <= 1)
                return false;

            PrecisedVertex h = new PrecisedVertex(-1e20d, v.Y, 0);
            PrecisedVertex w = new PrecisedVertex(v.X, v.Y, 0);
            int count = 0;

            if (!region[0].Equals(region[region.Count - 1]))
                region.Add(region[0]);

            for (int i = 1; i < region.Count; i++)
                if ((region[i - 1].X <= w.X || region[i].X <= w.X) && Between(w.Y, region[i - 1].Y, region[i].Y))
                {
                    PrecisedVertex cross = IntersectLineSegments2D(region[i - 1], region[i], h, w);
                    if (!cross.IsNull() && !region[i].Equals(cross) && (!w.Equals(cross) || (count & 1) == 0))
                        count++;
                }

            return (count & 1) != 0;
        }
        
        public static bool PointInRegion2DRevised(IVertexDouble v, IList<IVertexDouble> region)
        {
            if (region.Count <= 1)
                return false;

            PrecisedVertex h = new PrecisedVertex(-1e20d, v.Y, 0);
            PrecisedVertex w = new PrecisedVertex(v.X, v.Y, 0);
            int count = 0;

            if (!region[0].Equals(region[region.Count - 1]))
                region.Add(region[0]);

            for (int i = 0; i < region.Count; i++)
            {
                IVertexDouble previous = i == 0 ? region[region.Count - 1] : region[i - 1];

                if ((previous.X <= w.X || region[i].X <= w.X) && Between(w.Y, previous.Y, region[i].Y))
                {
                    PrecisedVertex cross = IntersectLineSegments2D(previous, region[i], h, w);
                    if (!cross.IsNull())
                    {
                        if(region[i].Equals(cross))
                        {
                            IVertexDouble next = i == region.Count - 1 ? region[0] : region[i + 1];
                            // попали на угол - пропускаем эту точку
                            if(IntersectLineSegments2D(previous, next, h, w).IsNull())
                                continue;
                        }
                        count++;
                    }
                }
            }

            return (count & 1) != 0;
        }

        public static bool IntersectsLineSegments2D(IVertexDouble Line1Point1, IVertexDouble Line1Point2,
                                                    IVertexDouble Line2Point1, IVertexDouble Line2Point2)
        {
            PrecisedVertex v = IntersectLines2D(Line1Point1, Line1Point2, Line2Point1, Line2Point2);
            return (!v.IsNull()) &&
                   (Line1Point1.X == Line1Point2.X || Between(v.X, Line1Point1.X, Line1Point2.X)) &&
                   (Line1Point1.Y == Line1Point2.Y || Between(v.Y, Line1Point1.Y, Line1Point2.Y)) &&
                   (Line2Point1.X == Line2Point2.X || Between(v.X, Line2Point1.X, Line2Point2.X)) &&
                   (Line2Point1.Y == Line2Point2.Y || Between(v.Y, Line2Point1.Y, Line2Point2.Y));
        }

        public static PrecisedVertex IntersectLineSegments2D(IVertexDouble Line1Point1, IVertexDouble Line1Point2,
                                                     IVertexDouble Line2Point1, IVertexDouble Line2Point2)
        {
            PrecisedVertex v = IntersectLines2D(Line1Point1, Line1Point2, Line2Point1, Line2Point2);
            return (!v.IsNull()) &&
                   (Line1Point1.X == Line1Point2.X || Between(v.X, Line1Point1.X, Line1Point2.X)) &&
                   (Line1Point1.Y == Line1Point2.Y || Between(v.Y, Line1Point1.Y, Line1Point2.Y)) &&
                   (Line2Point1.X == Line2Point2.X || Between(v.X, Line2Point1.X, Line2Point2.X)) &&
                   (Line2Point1.Y == Line2Point2.Y || Between(v.Y, Line2Point1.Y, Line2Point2.Y))
                       ? v
                       : PrecisedVertex.Null;
        }

        public static PrecisedVertex IntersectLines2D(IVertexDouble Line1Point1, IVertexDouble Line1Point2,
                                              IVertexDouble Line2Point1, IVertexDouble Line2Point2)
        {
            const float eps = 0.0000001f;
            double A1 = Line1Point1.Y - Line1Point2.Y;
            double B1 = Line1Point2.X - Line1Point1.X;
            double C1 = Line1Point2.Y*Line1Point1.X - Line1Point1.Y*Line1Point2.X;

            double A2 = Line2Point1.Y - Line2Point2.Y;
            double B2 = Line2Point2.X - Line2Point1.X;
            double C2 = Line2Point2.Y*Line2Point1.X - Line2Point1.Y*Line2Point2.X;

            double Det = A1*B2 - B1*A2;
            if (Math.Abs(Det) <= eps)
                return PrecisedVertex.Null;

            return new PrecisedVertex((B1*C2 - C1*B2)/Det, (C1*A2 - A1*C2)/Det, 0);
        }

        #endregion

        #region Другие функции


        public static double RadiansToDegrees(double radians)
        {
            return (radians / Math.PI) * 180;
        }

        public static double DegreesToRadians(double degrees)
        {
            return (degrees / 180) * Math.PI;
        }


        /// <summary>
        /// Изменить double.MaxValue на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = double.MaxValue, value в противном случае </returns>
        public static double ChangeDoubleMaxValueToValue(double value)
        {
            return ChangeDoubleMaxValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить double.MinValue на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = double.MinValue, value в противном случае </returns>
        public static double ChangeDoubleMinValueToValue(double value)
        {
            return ChangeDoubleMinValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить double.NaN на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = double.NaN, value в противном случае </returns>
        public static double ChangeDoubleNanValueToValue(double value)
        {
            return ChangeDoubleNanValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить double.MaxValue на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = double.MaxValue, value в противном случае </returns>
        public static double ChangeDoubleMaxValueToValue(double value, double newValue)
        {
            return value == double.MaxValue ? newValue : value;
        }
        /// <summary>
        /// Изменить double.MinValue на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = double.MinValue, value в противном случае </returns>
        public static double ChangeDoubleMinValueToValue(double value, double newValue)
        {
            return value == double.MinValue ? newValue : value;
        }
        /// <summary>
        /// Изменить double.NaN на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = double.NaN, value в противном случае </returns>
        public static double ChangeDoubleNanValueToValue(double value, double newValue)
        {
            return double.IsNaN(value) ? newValue : value;
        }

        /// <summary>
        /// Изменить float.MaxValue на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = float.MaxValue, value в противном случае </returns>
        public static float ChangeFloatMaxValueToValue(float value)
        {
            return ChangeFloatMaxValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить float.MinValue на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = float.MinValue, value в противном случае </returns>
        public static float ChangeFloatMinValueToValue(float value)
        {
            return ChangeFloatMinValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить float.NaN на ноль
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <returns>0, если значине = float.NaN, value в противном случае </returns>
        public static float ChangeFloatNanValueToValue(float value)
        {
            return ChangeFloatNanValueToValue(value, 0);
        }
        /// <summary>
        /// Изменить float.MaxValue на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = float.MaxValue, value в противном случае </returns>
        public static float ChangeFloatMaxValueToValue(float value, float newValue)
        {
            return value == float.MaxValue ? newValue : value;
        }
        /// <summary>
        /// Изменить float.MinValue на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = float.MinValue, value в противном случае </returns>
        public static float ChangeFloatMinValueToValue(float value, float newValue)
        {
            return value == float.MinValue ? newValue : value;
        }
        /// <summary>
        /// Изменить float.NaN на другое значение
        /// </summary>
        /// <param name="value">Проверяемое значение</param>
        /// <param name="newValue">Новое значение</param>
        /// <returns>newValue если значине = float.NaN, value в противном случае </returns>
        public static float ChangeFloatNanValueToValue(float value, float newValue)
        {
            return float.IsNaN(value) ? newValue : value;
        }
        /// <summary>
        /// Поворот на угол
        /// </summary>
        /// <param name="x">X,Y – координаты точки относительно которой происходит поворот</param>
        /// <param name="y">Y</param>
        /// <param name="alfa">Alfa –угол поворота по часовой стрелке</param>
        /// <param name="A">A – вертекс,который мы поворачиваем</param>
        public static Vertex RotationXYZ(double x, double y, double alfa, Vertex A)
        {
            double xx = x + (A.X - x) * Math.Cos(MathUtils.PI_180D * alfa) + (A.Y - y) * Math.Sin(MathUtils.PI_180D * alfa);
            double yy = y - (A.X - x) * Math.Sin(MathUtils.PI_180D * alfa) + (A.Y - y) * Math.Cos(MathUtils.PI_180D * alfa);
            A.Move((float)(xx), (float)(yy), A.Z);
            return A;
        }

        public static void RotateSurface(BaseSurface Surf, double alfa)
        {
            for (int i = 0; i < Surf.Nx; i++)
                for (int j = 0; j < Surf.Ny; j++)
                    if (Surf[0, 0].X != Surf.NullValue
                        && Surf[0, 0].Y != Surf.NullValue
                        && Surf[i, j].X != Surf.NullValue
                        && Surf[i, j].Y != Surf.NullValue
                        && Surf[i, j].Z != Surf.NullValue)
                    {
                        Surf[i, j] = RotationXYZ(Surf[0, 0].X, Surf[0, 0].Y, alfa, Surf[i, j]);
                    }
            Surf.RecalcMinMaxValue();
        }




        #endregion

        #region Площади

        public static double Striangle(Vertex A, Vertex B, Vertex C)
        {
            double x = (A.Y - B.Y) * (C.Z - B.Z) - (A.Z - B.Z) * (C.Y - B.Y);
            double y = (A.Z - B.Z) * (C.X - B.X) - (A.X - B.X) * (C.Z - B.Z);
            double z = (A.X - B.X) * (C.Y - B.Y) - (A.Y - B.Y) * (C.X - B.X);
            return Math.Sqrt(x * x + y * y + z * z) / 2;
        }
        public static double S2Dtriangle(Vertex A, Vertex B, Vertex C)
        {

            double z = (A.X - B.X) * (C.Y - B.Y) - (A.Y - B.Y) * (C.X - B.X);
            return Math.Abs(z) / 2;
        }

        private static double SofCellIJ(BaseSurface Surf, int i, int j)
        {
            double S = 0;
            if (Surf[i, j].Z != Surf.NullValue && Surf[i + 1, j + 1].Z != Surf.NullValue)
            {
                if (Surf[i + 1, j].Z != Surf.NullValue)
                    S += Striangle(Surf[i, j], Surf[i + 1, j], Surf[i + 1, j + 1]);
                if (Surf[i, j + 1].Z != Surf.NullValue)
                    S += Striangle(Surf[i, j], Surf[i, j + 1], Surf[i + 1, j + 1]);
                //continue;
            }else 
            if (Surf[i + 1, j].Z != Surf.NullValue && Surf[i, j + 1].Z != Surf.NullValue)
            {
                if (Surf[i, j].Z != Surf.NullValue)
                    S += Striangle(Surf[i, j + 1], Surf[i, j], Surf[i + 1, j]);
                if (Surf[i + 1, j + 1].Z != Surf.NullValue)
                    S += Striangle(Surf[i, j + 1], Surf[i + 1, j + 1], Surf[i + 1, j]);
                //continue;
            }
            return S;

        }
        private static double S2DofCellIJ(BaseSurface Surf, int i, int j)
        {
            double S = 0;
            if (Surf[i, j].Z != Surf.NullValue && Surf[i + 1, j + 1].Z != Surf.NullValue)
            {
                if (Surf[i + 1, j].Z != Surf.NullValue)
                    S += S2Dtriangle(Surf[i, j], Surf[i + 1, j], Surf[i + 1, j + 1]);
                if (Surf[i, j + 1].Z != Surf.NullValue)
                    S += S2Dtriangle(Surf[i, j], Surf[i, j + 1], Surf[i + 1, j + 1]);
                //continue;
            }
            else
                if (Surf[i + 1, j].Z != Surf.NullValue && Surf[i, j + 1].Z != Surf.NullValue)
                {
                    if (Surf[i, j].Z != Surf.NullValue)
                        S += S2Dtriangle(Surf[i, j + 1], Surf[i, j], Surf[i + 1, j]);
                    if (Surf[i + 1, j + 1].Z != Surf.NullValue)
                        S += S2Dtriangle(Surf[i, j + 1], Surf[i + 1, j + 1], Surf[i + 1, j]);
                    //continue;
                }
            return S;

        }
        public static double SSurface(BaseSurface Surf)
        {
            double S = 0;
            for (int i = 0; i < Surf.Nx - 1; i++)
                for (int j = 0; j < Surf.Ny - 1; j++)
                {
                    S += SofCellIJ(Surf, i, j); 
                }
            return S;
        }
        public static double S2DSurface(BaseSurface Surf)
        {
            double S = 0;
            for (int i = 0; i < Surf.Nx - 1; i++)
                for (int j = 0; j < Surf.Ny - 1; j++)
                {
                    S += S2DofCellIJ(Surf, i, j);
                }
            return S;
        }

        public static double SSurfaceInsidePolygon(BaseSurface Surf, List<TimeZYX.Model.Geology.Interpolation.PolygonInfo> polyInfo)
        {
            Geology.Interpolation.Color[,] active = Geology.Interpolation.Utils.ArrayToConsider(Surf, polyInfo);
            double S = 0; //Площадь

            //double V = 0;//Объем
            for (int i = 0; i < Surf.Nx - 1; i++)
                for (int j = 0; j < Surf.Ny - 1; j++)
                {
                    if (active[i, j].col % 2 == 1)
                    {
                        S += SofCellIJ(Surf, i, j);

                        //V += GetCellVolume(Surf, i, j, 0);
                    }
                }
            return S;
        }


        #endregion

        #region Объемы
    
        public static double VofPiramid(Vertex O, Vertex A, Vertex B, Vertex C, Vertex D)
        {
            double Sf = Striangle(A, B, C);
            double h = -1;
            if (Sf != 0) h = GetH(A, B, C, O);
            Sf += Striangle(C, D, A);
            if (Sf == 0) return 0;
            if (h < 0) h = GetH(C, D, A, O);
            return Sf * h / 3;
        }

        public static double VofPrisma(Vertex A, Vertex B, Vertex C, double h)
        {
            double Sf = Striangle(A, B, C);
            return Sf * h;
        }

        /// <summary>
        /// Объем пов-ти
        /// </summary>
        /// <param name="Surf"></param>
        /// <returns></returns>
        public static double VofSurface(BaseSurface Surf)
        {
            return VofSurfaceByLevel(Surf,0);
        }

        /// <summary>
        /// Объем пов-ти по уровню
        /// </summary>
        /// <param name="Surf"></param>
        /// <param name="Z0">Уровень</param>
        /// <returns></returns>
        public static double VofSurfaceByLevel(BaseSurface Surf, double Z0)
        {
            double V = 0;
            for (int i = 0; i < Surf.Nx - 1; i++)
                for (int j = 0; j < Surf.Ny - 1; j++)
                {
                     V += GetCellVolume(Surf, i, j, Z0);                
                }
            return Math.Round(V, 4); ;
        }

        public static double VofSurface(BaseSurface Surf1, BaseSurface Surf2)
        {           
            Vertex[,] newsurf = new Vertex[Surf1.Nx,Surf1.Ny];
            for (int i = 0; i < Surf1.Nx; i++)
                for (int j = 0; j < Surf1.Ny; j++)
                {
                    float z =(float) TimeZYX.Model.Geology.Interpolation.Utils.PointsPulling(Surf1[i, j].X, Surf1[i, j].Y, Surf2, Surf1.NullValue);
                    if (z != Surf1.NullValue) newsurf[i, j] = new Vertex(Surf1[i, j].X, Surf1[i, j].Y, Surf1[i, j].Z - z);
                    else newsurf[i, j] = new Vertex(Surf1[i, j].X, Surf1[i, j].Y, Surf1.NullValue);
                }
            BaseSurface surf = new Surface(Surf1.Project, "", Surf1.NullValue, newsurf);
            surf.RecalcMinMaxValue();
            return VofSurface(surf);
        }

        public static double VofSurfaceInsidePolygon(BaseSurface Surf, List<TimeZYX.Model.Geology.Interpolation.PolygonInfo> polyInfo)
        {
            Geology.Interpolation.Color[,] active = Geology.Interpolation.Utils.ArrayToConsider(Surf, polyInfo);
            //double S = 0; //Площадь

            double V = 0;//Объем
            for (int i = 0; i < Surf.Nx - 1; i++)
                for (int j = 0; j < Surf.Ny - 1; j++)
                {
                    if (active[i, j].col % 2 == 1)
                    {
                        //S += SofCellIJ(Surf, i, j);

                        V += GetCellVolume(Surf, i, j, 0);
                    }
                }
            return V;
        }
        static double VofTrianglePrizme(Vertex A0,Vertex B0, Vertex C0,double Z0)
        {   double nx,ny,nz,Jac,d,dV=0,u,v;
              Vertex m, l,A,C,B;
              if (A0.Z >= B0.Z)
                    {
                        if (A0.Z >= C0.Z)
                        {
                            B = A0;
                            A = B0;
                            C = C0;
                        }
                        else
                        {
                            B = C0;
                            A = A0;
                            C = B0;
                        }
                    }
              else if (B0.Z >= C0.Z)
                    {
                            B = B0;
                            A = A0;
                            C = C0;
                    }else 
                    {
                         B = C0;
                            A = A0;
                            C = B0;
                    }
         
                    if (B.Z - Z0 > 0)
                    {
                        TimeZYX.Model.Utils.WellTracerSDB.PlaneCoeff(1, A, B, C, out nx, out ny, out nz, out d);
                        if (nz != 0&&nz==nz)
                        {
                           
                            nx /= -nz;
                            ny /= -nz;
                            d /= -nz;
                            m = A - B;
                            l = C - B;
                            u =  (Z0 - B.Z) / m.Z;
                            if (m.Z == 0) u = 1;
                            else
                                if (u < 0) u = 0; else if (u > 1) u = 1;

                            v = (Z0 - B.Z) / l.Z;
                            if (l.Z == 0)  v = 1;
                            else
                                if (v < 0) v = 0; else if (v > 1) v = 1;
                            
                            Jac = Math.Abs(m.X * l.Y - m.Y * l.X) / 2;
                            dV = Jac * ((B.X + (m.X*u*u + m.Y*v*v) * 0.5) * nx + (B.Y + 0.5 * (l.X*u*u + l.Y*v*v)) * ny + d-Z0)*u*v;
                            if (dV < 0) dV = 0 ;
                           
                        }
                    }
                    return dV;

        }
        static double GetCellVolume(BaseSurface Surf, int i, int j, double Z0)
        {
            double V = 0,dV;
          
            if (Surf[i, j].Z != Surf.NullValue && Surf[i + 1, j + 1].Z != Surf.NullValue)
            {
                if (Surf[i + 1, j].Z != Surf.NullValue)
                {
                    dV = VofTrianglePrizme(Surf[i, j], Surf[i+1, j], Surf[i+1, j+1], Z0);
                    V += dV;
                  
                }
                if (Surf[i, j + 1].Z != Surf.NullValue)
                {
                    dV = VofTrianglePrizme(Surf[i, j], Surf[i, j+1], Surf[i + 1, j + 1], Z0);
                    V += dV;
                }
                //continue;
            }else
                if (Surf[i + 1, j].Z != Surf.NullValue && Surf[i, j + 1].Z != Surf.NullValue)
                {
                    if (Surf[i, j].Z != Surf.NullValue)
                    {
                        dV = VofTrianglePrizme(Surf[i+1, j], Surf[i, j], Surf[i , j + 1], Z0);
                         V += dV;
                    }
                    if (Surf[i + 1, j + 1].Z != Surf.NullValue)
                    {
                        dV = VofTrianglePrizme(Surf[i + 1, j], Surf[i, j], Surf[i, j + 1], Z0);
                        V += dV;
                    }
                }

          
            return V;

        }

        static double GetH(Vertex O, Vertex A, Vertex B, Vertex C)
        {
            double x = (A.Y - B.Y) * (C.Z - B.Z) - (A.Z - B.Z) * (C.Y - B.Y);
            double y = (A.Z - B.Z) * (C.X - B.X) - (A.X - B.X) * (C.Z - B.Z);
            double z = (A.X - B.X) * (C.Y - B.Y) - (A.Y - B.Y) * (C.X- B.X);
            double s = Math.Sqrt(x * x + y * y + z * z);

            if (s == 0)
                return 0;

            x /= s;
            y /= s;
            z /= s;
            return (Math.Abs(x * (B.X - O.X) + y * (B.Y - O.Y) + z * (B.Z - O.Z)));
            //double x = (A.Y - B.Y) * (C.Z - B.Z) - (A.Z - B.Z) * (C.Y - B.Y);
            //double y = (A.Z - B.Z) * (C.X - B.X) - (A.X - B.X) * (C.Z - B.Z); double z = (A.X - B.X) * (C.Y - B.Y) - (A.Y - B.Y) * (C.Y - B.Y);
            //double s = Math.Sqrt(x * x + y * y + z * z);
            //x /= s;
            //y /= s;
            //z /= s;
            //return (Math.Abs(x * (B.X - O.X) + y * (B.Y - O.Y) + z * (B.Z - O.Z)));
        }


        #endregion

        /// <summary>
        /// Проверка, с какой стороны от прямой (p1, p2) оказалась точка p
        /// </summary>
        /// <param name="p1">Точка, задающая прямую</param>
        /// <param name="p2">Точка, задающая прямую</param>
        /// <param name="p">Проверяемая точка</param>
        /// <returns>0, если точка на прямой, -1 и 1 в зависимости от стороны, с которой точка оказалась от прямой. Точки с одной стороны дадут один и тот же результат</returns>
        private static int GetLineSideXY(IVertexFloat p1, IVertexFloat p2, IVertexFloat p)
        {
            return Math.Sign((p2.X - p1.X) * (p.Y - p1.Y) - (p2.Y - p1.Y) * (p.X - p1.X));
        }

        /// <summary>
        /// Проверка, образуют ли данные 4 точки выпуклый четырехугольник
        /// </summary>
        /// <param name="p0">Точка 1</param>
        /// <param name="p1">Точка 2</param>
        /// <param name="p2">Точка 3</param>
        /// <param name="p3">Точка 4</param>
        /// <returns>true, если точки представляют собой собственную выпуклую оболочку</returns>
        public static bool IsConvexTetragon(IVertexFloat p0, IVertexFloat p1, IVertexFloat p2, IVertexFloat p3)
        {
            bool isConvex = true;

            isConvex &= (GetLineSideXY(p0, p1, p2) * GetLineSideXY(p0, p1, p3) > 0);
            isConvex &= (GetLineSideXY(p1, p2, p3) * GetLineSideXY(p1, p2, p0) > 0);
            isConvex &= (GetLineSideXY(p2, p3, p0) * GetLineSideXY(p2, p3, p1) > 0);
            isConvex &= (GetLineSideXY(p3, p0, p1) * GetLineSideXY(p3, p0, p2) > 0);

            return isConvex;
        }

        public static double Bfunc(double tj, double tj1, double tj2, double tj3, double t)
        {
            double b = 0;
            if (t < tj || t > tj3) return 0;
            if (t < tj1)
            {
                b = (t - tj) * (t - tj) / ((tj2 - tj) * (tj1 - tj));
                return (b);
            }
            if (t < tj2)
            {
                b = (t - tj) * (tj2 - t) / ((tj2 - tj) * (tj2 - tj1)) + (tj3 - t) * (t - tj1) / ((tj3 - tj1) * (tj2 - tj1));
                return (b);
            }
            b = (tj3 - t) * (tj3 - t) / ((tj3 - tj1) * (tj3 - tj2));

            return (b);
        }

        public static void GetSpline(int TP, double[] Li, PrecisedVertex[] data, out double[] px, out double[] py, out double[] pz)
        {
            px = new double[1]; py = new double[1]; pz = new double[1];
            int M = -1;
            int count = 1;
            for (int i = 1; i < Li.Length; i++)
            {
                if (Math.Abs(Li[i] - Li[i - 1]) > 0.00001) count++;
            }
            if (count == 1) {  return; }
         
            if (TP == 0)
            {
                M = 2;
                double[,] A = new double[M, M];
                double[] c = new double[M];

                double[] b = new double[M];
                px = new double[M];
                py = new double[M];
                pz = new double[M];

                for (int i = 0; i < data.Length; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        c[j] = Math.Pow(Li[i] / Li[Li.Length - 1], j);

                    }
                    for (int j = 0; j < M; j++)
                    {

                        b[j] += c[j] * data[i].Z;
                        for (int k = 0; k < M; k++)
                        {
                            A[j, k] += c[j] * c[k];

                        }
                    }
                    px[0] += data[i].X;
                    py[0] += data[i].Y;
                }
                px[0] /= data.Length;
                py[0] /= data.Length;
                MathFunctions.GaussSolve(ref A, ref b, ref pz, M);
                return;
            }
             M = 2;
         
         
             if (TP == 2 && count > 2) M = 3;
             if (TP == 3 && count > 4) M = 5;
            double[,] Ax = new double[M, M];
            double[,] Ay = new double[M, M];
            double[,] Az = new double[M, M];
            double[] coef = new double[M];
            double[] bx = new double[M];
            double[] by = new double[M];
            double[] bz = new double[M];
            px = new double[M];
            py = new double[M];
            pz = new double[M];

            for (int i = 0; i < data.Length; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    coef[j] = Math.Pow(Li[i] / Li[Li.Length - 1], j);

                }
                for (int j = 0; j < M; j++)
                {
                    bx[j] += coef[j] * data[i].X;
                    by[j] += coef[j] * data[i].Y;
                    bz[j] += coef[j] * data[i].Z;
                    for (int k = 0; k < M; k++)
                    {
                        Ax[j, k] += coef[j] * coef[k];

                    }
                }

            }
            Array.Copy(Ax, Ay, Ax.Length);
            Array.Copy(Ax, Az, Ax.Length);
            MathFunctions.GaussSolve(ref Ax, ref bx, ref px, M);
            MathFunctions.GaussSolve(ref Ay, ref by, ref py, M);
            MathFunctions.GaussSolve(ref Az, ref bz, ref pz, 2);




        }

        public static void GaussSolve(ref double[,] a, ref double[] b, ref double[] x, int neq)
        {
            int i, j, k, kp1;
            int[] ipvt = new int[neq];
            double[] y = new double[neq];

            if (neq == 1)
            {
                x[0] = b[0] / a[0, 0];
                return;
            }

            // change consequence on matrix lines
            double help = 0;


            // Decompose
            //ipvt[neq - 1] = 1;

            double t;
            /*  gaussian elimination  */
            for (k = 1; k < neq; ++k)
            {
                kp1 = k + 1;
                bool flag = false;
                if (a[k - 1, k - 1] == 0)
                {

                    for (j = k + 1; j <= neq; ++j)
                    {

                        if (a[j - 1, k - 1] != 0)
                        {
                            for (i = 0; i < neq; ++i)
                            {
                                y[i] = a[k - 1, i];
                                a[k - 1, i] = a[j - 1, i];
                                a[j - 1, i] = y[i];
                            }
                            help = b[k - 1];
                            b[k - 1] = b[j - 1];
                            b[j - 1] = help;
                            flag = true;
                            break;
                        }
                    }
                    if (!flag) { throw new Exception(Properties.Resources.Const_115); }

                }
                for (j = kp1; j <= neq; ++j)
                {
                    if (a[j - 1, k - 1] == 0.0) continue;
                    t = -a[k - 1, k - 1] / a[j - 1, k - 1];

                    for (i = k + 1; i <= neq; ++i)
                    {
                        if (a[j - 1, i - 1] == a[j - 1, k - 1])
                        {
                            a[j - 1, i - 1] = (a[k - 1, i - 1] - a[k - 1, k - 1]);
                            continue;
                        }
                        a[j - 1, i - 1] *= t;
                        a[j - 1, i - 1] = a[k - 1, i - 1] + a[j - 1, i - 1];

                    }
                    a[j - 1, k - 1] = 0;
                    b[j - 1] *= t;
                    b[j - 1] = b[k - 1] + b[j - 1];
                }
            }

            /* back substitution */
            for (k = neq; k > 0; --k)
            {
                help = 0.0f;
                for (i = k + 1; i <= neq; ++i) help += b[i - 1] * a[k - 1, i - 1];//re
                x[k - 1] = (b[k - 1] - help) / a[k - 1, k - 1];
            }
        }

        public static PrecisedVertex GetXYZ(double tt, double a, double b, double[] px, double[] py, double[] pz)
        {
            PrecisedVertex V;
            double t = tt * a + b;
            double x, y, z;
            x = y = z = 0;
            if (t < 0)
            {
                x = px[1] * t + px[0];
                y = py[1] * t + py[0];
                z = pz[1] * t + pz[0];

                V = new PrecisedVertex(x, y, z);
                return (V);
            }
            if (t > 1)
            {

                x = px[0];
                y = py[0];
                z = pz[0];
                for (int j = 1; j < px.Length; j++)
                {
                    x += px[j] * ((t - 1) * j + 1);
                    y += py[j] * ((t - 1) * j + 1);
                    z += pz[j] * ((t - 1) * j + 1);
                }
                V = new PrecisedVertex(x,y, z);
                return (V);
            }
            for (int j = 0; j < px.Length; j++)
            {
                x += px[j] * Math.Pow(t, j);
                y += py[j] * Math.Pow(t, j);
                z += pz[j] * Math.Pow(t, j);
            }
            V = new PrecisedVertex(x, y, z);
            return (V);
        }

        public static void GetTNormCoeff( double Zbottom, double Ztop, double[] px, double[] py, double[] pz, out double a, out double b)
        {
            a = 1;
            b = 0;
        
          
                List<double> Roots = new List<double>(2);
                List<double> buf = new List<double>(4);
                double t, az0, azN, z0, zN;
                z0 = pz[0];
                zN = pz[0];
                az0 = -pz[1];
                azN = 0;
                for (int j = 1; j < pz.Length; j++)
                {
                    azN += j * pz[j];
                    zN += pz[j];
                }
                if (az0 == 0) Roots.Add(0);
                else
                {
                    t = (Zbottom - z0) / az0;
                    if (t >= 0) Roots.Add(-t);
                    t = (Ztop - z0) / az0;

                    if (t >= 0) Roots.Add(-t);
                }
                if (azN == 0) Roots.Add(1);
                else
                {
                    t = (Zbottom - zN) / azN;
                    if (t >= 0) Roots.Add(t + 1);
                    t = (Ztop - zN) / azN;
                    if (t >= 0) Roots.Add(t + 1);
                }
            if (Roots.Count<2)
                switch (pz.Length)
                {
                    case 2:
                        t = (Zbottom - z0) / (zN - z0);
                        if (t >= 0 && t <= 1) Roots.Add(t);
                        t = (Ztop - z0) / (zN - z0);
                        if (t >= 0 && t <= 1) Roots.Add(t);
                        break;
                    case 3:
                        double des = pz[1] * pz[1] - 4 * pz[2] * (pz[0] - Ztop);
                        if (des > 0)
                        {
                            t = (-pz[1] + Math.Sqrt(des)) / 2 / pz[2];
                            if (t > 0 && t < 1) Roots.Add(t);
                            t = (-pz[1] - Math.Sqrt(des)) / 2 / pz[2];
                            if (t > 0 && t < 1) Roots.Add(t);
                        }
                        des = pz[1] * pz[1] - 4 * pz[2] * (pz[0] - Zbottom);
                        if (des > 0)
                        {
                            t = (-pz[1] + Math.Sqrt(des)) / 2 / pz[2];
                            if (t >= 0 && t <= 1) Roots.Add(t);
                            t = (-pz[1] - Math.Sqrt(des)) / 2 / pz[2];
                            if (t >= 0 && t <= 1) Roots.Add(t);
                        }
                        break;
                    case 5:

                        buf = QuartEquationSolution(pz[4], pz[3], pz[2], pz[1], pz[0] - Zbottom);
                        for (int i = 0; i < buf.Count; i++)
                            if (buf[i] > 0 && buf[i] < 1)
                                Roots.Add(buf[i]);
                        buf = QuartEquationSolution(pz[4], pz[3], pz[2], pz[1], pz[0] - Ztop);
                        for (int i = 0; i < buf.Count; i++)
                            if (buf[i] > 0 && buf[i] < 1)
                                Roots.Add(buf[i]);
                        break;
                }
                if (Roots.Count < 2) return;
                double[] forSort = new double[Roots.Count];
                double[] roots = Roots.ToArray();
                for (int i = 0; i < forSort.Length; i++)
                    forSort[i] = Math.Abs(Roots[i] - 0.5);
                Array.Sort(forSort, roots);
                if (roots[0] > roots[1])
                {
                    b = roots[1];
                    a = roots[0];
                }
                else
                {
                    b = roots[0];
                    a = roots[1];
                }
                return;
        
        }

        public static List<double> QuartEquationSolution(double a, double b, double c, double d, double e)
        {

            List<double> x = new List<double>(4);
            double alfa = -3 * b * b / (8 * a * a) + c / a;
            double betta = b * b * b / (8 * a * a * a) - b * c / (2 * a * a) + d / a;
            double gamma = -3 * b * b * b * b / (256 * a * a * a * a) + c * b * b / (16 * a * a * a) - b * d / (4 * a * a) + e / a;
            double buf;
            Complex X;
            if (betta == 0)
            {
                buf = -b / (4 * a);
                X = -buf + MathUtils.CSqrt((-alfa + MathUtils.CSqrt(alfa * alfa - 4 * gamma))) / 2;
                if (X.y == 0) x.Add(X.x);
                X = -buf - MathUtils.CSqrt((-alfa + MathUtils.CSqrt(alfa * alfa - 4 * gamma))) / 2;
                if (X.y == 0) x.Add(X.x);
                X = -buf + MathUtils.CSqrt((-alfa - MathUtils.CSqrt(alfa * alfa - 4 * gamma))) / 2;
                if (X.y == 0) x.Add(X.x);
                X = -buf - MathUtils.CSqrt((-alfa - MathUtils.CSqrt(alfa * alfa - 4 * gamma))) / 2;
                if (X.y == 0) x.Add(X.x);
                return x;

            }
            Complex P = -alfa * alfa / 12 - gamma;
            Complex Q = -alfa * alfa * alfa / 108 + alfa * gamma / 3 - betta * betta / 8;


            Complex R = MathUtils.CSqrtn(-Q / 2 - MathUtils.CSqrt(Q * Q / 4 + P * P * P / 27), 3);
            // Complex R2 = -new Complex(Math.Pow(-(-Q / 2 - MathUtils.CSqrt(Q * Q / 4 + P * P * P / 27)).x, 1.0/3),0);
            Complex Y;
            //    Complex R = R1;
            if (R.x == 0 && R.y == 0)
                Y = -5 * alfa / 6 + R - Math.Pow(Q.x, 1.0 / 3);
            else
                Y = -5 * alfa / 6 + R - P / (3 * R);
            buf = -b / (4 * a);
            Complex W = MathUtils.CSqrt(alfa + 2 * Y);
            X = buf + 0.5 * (W + MathUtils.CSqrt(-(3 * alfa + 2 * Y + 2 * betta / W)));
            if (Math.Abs(X.y) < 0.000001) x.Add(X.x);
            X = buf + 0.5 * (W - MathUtils.CSqrt(-(3 * alfa + 2 * Y + 2 * betta / W)));
            if (Math.Abs(X.y) < 0.000001) x.Add(X.x);
            X = buf + 0.5 * (-W + MathUtils.CSqrt(-(3 * alfa + 2 * Y - 2 * betta / W)));
            if (Math.Abs(X.y) < 0.000001) x.Add(X.x);
            X = buf + 0.5 * (-W - MathUtils.CSqrt(-(3 * alfa + 2 * Y - 2 * betta / W)));
            if (Math.Abs(X.y) < 0.000001) x.Add(X.x);
            return x;
        }

        public static double[] GetLengtSpline(PrecisedVertex[] data)
        {
            double[] Li = new double[data.Length];
            Li[0] = 0;
            for (int i = 1; i < data.Length; i++)
                Li[i] = Li[i - 1] + Math.Sqrt((data[i].X - data[i - 1].X) * (data[i].X - data[i - 1].X) + (data[i].Y - data[i - 1].Y) * (data[i].Y - data[i - 1].Y) + (data[i].Z - data[i - 1].Z) * (data[i].Z - data[i - 1].Z));
            return Li;
        }



        #region Безье
        /// <summary>
        /// Получение контрольных точек для построения кривой Безье
        /// </summary>
        /// <param name="data">Исходные координаты</param>
        /// <param name="P1">Контрольные точки</param>
        /// <param name="P2">Контрольные точки</param>
        public static void GetControlPoints(PrecisedVertex[] data, out PrecisedVertex[] P1, out PrecisedVertex[] P2)
        {
            PrecisedVertex r1, r2, Pbuf;
            P1 = new PrecisedVertex[data.Length - 1];
            P2 = new PrecisedVertex[data.Length - 1];
            double DerivateParam =0.2;
            int k = 0;
            double l = 0,lZ=0;
            double r1module = 0;
            double r2module = 0;
            double rmodule = 0;
            if (data.Length==2)
            {
                P1[0] = new PrecisedVertex(data[0].X + 0.25 * (data[1].X - data[0].X), data[0].Y+ 0.25 * (data[1].Y - data[0].Y), data[0].Z + 0.25 * (data[1].Z - data[0].Z));
                P2[0] = new PrecisedVertex(data[0].X + 0.75 * (data[1].X - data[0].X), data[0].Y + 0.75 * (data[1].Y - data[0].Y), data[0].Z + 0.75 * (data[1].Z - data[0].Z));
                return;
            }
            r1 = new PrecisedVertex(data[k].X - data[k + 1].X, data[k].Y - data[k + 1].Y, data[k].Z - data[k + 1].Z);
            r1module = Math.Sqrt(r1.X * r1.X + r1.Y * r1.Y + r1.Z * r1.Z);
            r2 = new PrecisedVertex(data[k + 2].X - data[k + 1].X, data[k + 2].Y - data[k + 1].Y, data[k + 2].Z - data[k + 1].Z);
            r2module = Math.Sqrt(r2.X * r2.X + r2.Y * r2.Y + r2.Z * r2.Z);
            Pbuf = new PrecisedVertex(r1.X / r1module + r2.X / r2module, r1.Y / r1module + r2.Y / r2module, r1.Z / r1module + r2.Z/r2module);
          
            Pbuf = Pbuf.CrossProduct(r1.CrossProduct(r2));
            rmodule=Math.Sqrt(Pbuf.X * Pbuf.X + Pbuf.Y * Pbuf.Y + Pbuf.Z * Pbuf.Z);
            lZ = 0;
            if (Math.Abs(rmodule) < 0.0001) 
                l = 0;
            else
            {
                Pbuf = new PrecisedVertex(Pbuf.X / rmodule, Pbuf.Y / rmodule, Pbuf.Z / rmodule);
                l = DerivateParam * (Pbuf.X * r1.X + Pbuf.Y * r1.Y + Pbuf.Z * r1.Z);
                if (Math.Abs(Pbuf.Z) > 0.00001)
                {
                    lZ = ((r1.Z) / Pbuf.Z);
                    if (lZ < 0) 
                        lZ = 0; 
                    else if (lZ >= l) 
                        lZ = l;
                }
                else 
                    lZ = l;
                
            }
            P2[k] = new PrecisedVertex(Pbuf.X * l + data[k + 1].X, Pbuf.Y * l + data[k + 1].Y, Pbuf.Z * lZ + data[k + 1].Z);
            P1[k] = new PrecisedVertex(0.5 * (P2[k].X + data[k].X), 0.5 * (P2[k].Y + data[k].Y), 0.5 * (P2[k].Z + data[k].Z));
            
            
            
            
            for (k = 1; k < P1.Length - 1; k++)
            {
                r1 = new PrecisedVertex(data[k].X - data[k + 1].X, data[k].Y - data[k + 1].Y, data[k].Z - data[k + 1].Z);
                r1module = Math.Sqrt(r1.X * r1.X + r1.Y * r1.Y + r1.Z * r1.Z);
                r2 = new PrecisedVertex(data[k + 2].X - data[k + 1].X, data[k + 2].Y - data[k + 1].Y, data[k + 2].Z - data[k + 1].Z);
                r2module = Math.Sqrt(r2.X * r2.X + r2.Y * r2.Y + r2.Z * r2.Z);
                rmodule = Math.Sqrt(Pbuf.X * Pbuf.X + Pbuf.Y * Pbuf.Y + Pbuf.Z * Pbuf.Z);
                lZ = 0;
                if (Math.Abs(rmodule) < 0.0001) 
                    l = 0;
                else
                {
                    Pbuf = new PrecisedVertex(Pbuf.X / rmodule, Pbuf.Y / rmodule, Pbuf.Z / rmodule);
                    l = DerivateParam * (Pbuf.X * r1.X + Pbuf.Y * r1.Y + Pbuf.Z * r1.Z);
                    if (Math.Abs(Pbuf.Z) > 0.00001)
                    {
                        lZ = ((r1.Z) / Pbuf.Z);
                        if (lZ < 0)
                            lZ = 0;
                        else if (lZ >= l) 
                            lZ = l;
                    }
                    else 
                        lZ = l;
                }
                P1[k] = new PrecisedVertex(-Pbuf.X*l + data[k].X, -Pbuf.Y*l + data[k].Y, -Pbuf.Z*lZ + data[k].Z);
                Pbuf = new PrecisedVertex(r1.X / r1module + r2.X / r2module, r1.Y / r1module + r2.Y / r2module, r1.Z / r1module + r2.Z / r2module);
            
                Pbuf = Pbuf.CrossProduct(r1.CrossProduct(r2));
                rmodule = Math.Sqrt(Pbuf.X * Pbuf.X + Pbuf.Y * Pbuf.Y + Pbuf.Z * Pbuf.Z);
                lZ = 0;
                if (Math.Abs(rmodule) < 0.0001) l = 0;
                else
                {
                    Pbuf = new PrecisedVertex(Pbuf.X / rmodule, Pbuf.Y / rmodule, Pbuf.Z / rmodule);
                    l = DerivateParam * (Pbuf.X * r1.X + Pbuf.Y * r1.Y + Pbuf.Z * r1.Z);
                    if (Math.Abs(Pbuf.Z) > 0.00001)
                    {
                        lZ = ((r1.Z) / Pbuf.Z);
                        if (lZ < 0) lZ = 0; else if (lZ >= l) lZ = l;
                    }
                    else lZ = l;
                }
                P2[k] = new PrecisedVertex(Pbuf.X * l + data[k + 1].X, Pbuf.Y * l + data[k + 1].Y, Pbuf.Z * lZ + data[k + 1].Z);

            }
            r1 = new PrecisedVertex(data[k].X - data[k + 1].X, data[k].Y - data[k + 1].Y, data[k].Z - data[k + 1].Z);
            rmodule = Math.Sqrt(Pbuf.X * Pbuf.X + Pbuf.Y * Pbuf.Y + Pbuf.Z * Pbuf.Z);
            lZ = 0;
            if (Math.Abs(rmodule) < 0.0001) l = 0;
            else
            {
                Pbuf = new PrecisedVertex(Pbuf.X / rmodule, Pbuf.Y / rmodule, Pbuf.Z / rmodule);
                l = DerivateParam * (Pbuf.X * r1.X + Pbuf.Y * r1.Y + Pbuf.Z * r1.Z);
                if (Math.Abs(Pbuf.Z) > 0.00001)
                {
                    lZ = ((r1.Z) / Pbuf.Z);
                    if (lZ < 0) lZ = 0; else if (lZ >= l) lZ = l;
                }
                else lZ = l;
            }
            P1[k] = new PrecisedVertex(-Pbuf.X * l + data[k].X, -Pbuf.Y * l + data[k].Y, -Pbuf.Z * lZ + data[k].Z);
            P2[k] = new PrecisedVertex(0.5 * (P1[k].X + data[data.Length - 1].X), 0.5 * (P1[k].Y + data[data.Length - 1].Y), 0.5 * (P1[k].Z + data[data.Length - 1].Z));
        }

        /// <summary>
        /// Построение кривой Безье
        /// </summary>
        /// <param name="P0">Первая точка</param>
        /// <param name="P1">Первая контрольная точка</param>
        /// <param name="P2">Вторая контрольная точка</param>
        /// <param name="P3">Вторая точка</param>
        /// <param name="t">Коэфициент Безье(шаг)</param>
        /// <returns>Новая точка кривой Безье</returns>
        public static PrecisedVertex BezierCurve(PrecisedVertex P0, PrecisedVertex P1, PrecisedVertex P2, PrecisedVertex P3, double t)
        {
            double x = P0.X * Math.Pow((1 - t), 3) + 3 * Math.Pow((1 - t), 2) * t * P1.X + 3 * (1 - t) * t * t * P2.X + t * t * t * P3.X;
            double y = P0.Y * Math.Pow((1 - t), 3) + 3 * Math.Pow((1 - t), 2) * t * P1.Y + 3 * (1 - t) * t * t * P2.Y + t * t * t * P3.Y;
            double z = P0.Z * Math.Pow((1 - t), 3) + 3 * Math.Pow((1 - t), 2) * t * P1.Z + 3 * (1 - t) * t * t * P2.Z + t * t * t * P3.Z;
            return (new PrecisedVertex((float)x, (float)y, (float)z));
        } 
        #endregion
        #region Массивы
        public static float GetArrayMin(float[] a)
        {
            int NCompare=a.Length/2;
            int Dind =1;
            int Dind_next = 2;
            int Poor = -1;
            int k,i;
            float t;
            int[] indexes = new int[a.Length];
                for (i = 0; i < a.Length; i++)
                {
                    indexes[i]=i;
                }
            while (NCompare>0)
            {
                Dind_next=Dind*2;
                for (i = 0; i < NCompare; i++)
                {
                    k=Dind_next*i;
                    if (a[indexes[k]] > a[indexes[k + 1]]) { indexes[k] = indexes[k + 1]; }
                }
                Dind = Dind_next;
                NCompare = a.Length / Dind;
            }
            return indexes[0];
        }
        
        #endregion
    }

    public struct Complex
    {
        public double x;
        public double y;

        public Complex(double _x)
        {
            x = _x;
            y = 0;
        }
        public Complex(double _x, double _y)
        {
            x = _x;
            y = _y;
        }
        public static implicit operator Complex(double _x)
        {
            return new Complex(_x);
        }
        public static bool operator ==(Complex lhs, Complex rhs)
        {
            return (lhs.x == rhs.x) & (lhs.y == rhs.y);
        }
        public static bool operator !=(Complex lhs, Complex rhs)
        {
            return (lhs.x != rhs.x) | (lhs.y != rhs.y);
        }
        public static Complex operator +(Complex lhs)
        {
            return lhs;
        }
        public static Complex operator -(Complex lhs)
        {
            return new Complex(-lhs.x, -lhs.y);
        }
        public static Complex operator +(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x + rhs.x, lhs.y + rhs.y);
        }
        public static Complex operator -(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x - rhs.x, lhs.y - rhs.y);
        }
        public static Complex operator *(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x * rhs.x - lhs.y * rhs.y, lhs.x * rhs.y + lhs.y * rhs.x);
        }
        public static Complex operator /(Complex lhs, Complex rhs)
        {
            Complex result;
            double e;
            double f;
            if (System.Math.Abs(rhs.y) < System.Math.Abs(rhs.x))
            {
                e = rhs.y / rhs.x;
                f = rhs.x + rhs.y * e;
                result.x = (lhs.x + lhs.y * e) / f;
                result.y = (lhs.y - lhs.x * e) / f;
            }
            else
            {
                e = rhs.x / rhs.y;
                f = rhs.y + rhs.x * e;
                result.x = (lhs.y + lhs.x * e) / f;
                result.y = (-lhs.x + lhs.y * e) / f;
            }
            return result;
        }

        public bool Equals(Complex other)
        {
            return other.x == x && other.y == y;
        }

        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (obj.GetType() != typeof(Complex)) return false;
            return Equals((Complex)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                return (x.GetHashCode() * 397) ^ y.GetHashCode();
            }
        }
    }

    public static class Fourier
    {
        /*************************************************************************
    Быстрое преобразование Фурье

    Алгоритм проводит быстрое преобразование Фурье вещественной
    функции, заданной n отсчетами на действительной оси.

    В зависимости от  переданных параметров, может выполняться
    как прямое, так и обратное преобразование.

    Входные параметры:
        tnn  -   Число значений функции. Должно  быть  степенью
                двойки. Алгоритм   не  проверяет  правильность
                переданного значения.
        a   -   array [0 .. nn-1] of Real
                Значения функции.
        InverseFFT
            -   направление преобразования.
                True, если обратное, False, если прямое.
                
    Выходные параметры:
        a   -   результат   преобразования.   Подробнее    см.
                описание на сайте.
    *************************************************************************/
        public static void RealFastFourierTransform(ref double[] a,
                                                    int tnn,
                                                    bool inversefft)
        {
            double twr;
            double twi;
            double twpr;
            double twpi;
            double twtemp;
            double ttheta;
            int i;
            int i1;
            int i2;
            int i3;
            int i4;
            double c1;
            double c2;
            double h1r;
            double h1i;
            double h2r;
            double h2i;
            double wrs;
            double wis;
            int ii;
            int m;
            int isign;
            double tempr;
            double tempi;

            if (tnn == 1)
            {
                return;
            }
            if (!inversefft)
            {
                ttheta = 2 * System.Math.PI / tnn;
                c1 = 0.5;
                c2 = -0.5;
            }
            else
            {
                ttheta = 2 * System.Math.PI / tnn;
                c1 = 0.5;
                c2 = 0.5;
                ttheta = -ttheta;
                twpr = -(2.0 * MathFunctions.Pow2(System.Math.Sin(0.5 * ttheta)));
                twpi = System.Math.Sin(ttheta);
                twr = 1.0 + twpr;
                twi = twpi;
                for (i = 2; i <= tnn / 4 + 1; i++)
                {
                    i1 = i + i - 2;
                    i2 = i1 + 1;
                    i3 = tnn + 1 - i2;
                    i4 = i3 + 1;
                    wrs = twr;
                    wis = twi;
                    h1r = c1 * (a[i1] + a[i3]);
                    h1i = c1 * (a[i2] - a[i4]);
                    h2r = -(c2 * (a[i2] + a[i4]));
                    h2i = c2 * (a[i1] - a[i3]);
                    a[i1] = h1r + wrs * h2r - wis * h2i;
                    a[i2] = h1i + wrs * h2i + wis * h2r;
                    a[i3] = h1r - wrs * h2r + wis * h2i;
                    a[i4] = -h1i + wrs * h2i + wis * h2r;
                    twtemp = twr;
                    twr = twr * twpr - twi * twpi + twr;
                    twi = twi * twpr + twtemp * twpi + twi;
                }
                h1r = a[0];
                a[0] = c1 * (h1r + a[1]);
                a[1] = c1 * (h1r - a[1]);
            }
            if (inversefft)
            {
                isign = -1;
            }
            else
            {
                isign = 1;
            }
            int n = tnn;
            int nn = tnn / 2;
            int j = 1;
            for (ii = 1; ii <= nn; ii++)
            {
                i = 2 * ii - 1;
                if (j > i)
                {
                    tempr = a[j - 1];
                    tempi = a[j];
                    a[j - 1] = a[i - 1];
                    a[j] = a[i];
                    a[i - 1] = tempr;
                    a[i] = tempi;
                }
                m = n / 2;
                while (m >= 2 & j > m)
                {
                    j = j - m;
                    m = m / 2;
                }
                j = j + m;
            }
            int mmax = 2;
            while (n > mmax)
            {
                int istep = 2 * mmax;
                double theta = 2 * System.Math.PI / (isign * mmax);
                double wpr = -(2.0 * MathFunctions.Pow2(System.Math.Sin(0.5 * theta)));
                double wpi = System.Math.Sin(theta);
                double wr = 1.0;
                double wi = 0.0;
                for (ii = 1; ii <= mmax / 2; ii++)
                {
                    m = 2 * ii - 1;
                    int jj;
                    for (jj = 0; jj <= (n - m) / istep; jj++)
                    {
                        i = m + jj * istep;
                        j = i + mmax;
                        tempr = wr * a[j - 1] - wi * a[j];
                        tempi = wr * a[j] + wi * a[j - 1];
                        a[j - 1] = a[i - 1] - tempr;
                        a[j] = a[i] - tempi;
                        a[i - 1] = a[i - 1] + tempr;
                        a[i] = a[i] + tempi;
                    }
                    double wtemp = wr;
                    wr = wr * wpr - wi * wpi + wr;
                    wi = wi * wpr + wtemp * wpi + wi;
                }
                mmax = istep;
            }
            if (inversefft)
            {
                for (i = 1; i <= 2 * nn; i++)
                {
                    a[i - 1] = a[i - 1] / nn;
                }
            }
            if (inversefft) return;
            twpr = -(2.0 * MathFunctions.Pow2(System.Math.Sin(0.5 * ttheta)));
            twpi = System.Math.Sin(ttheta);
            twr = 1.0 + twpr;
            twi = twpi;
            for (i = 2; i <= tnn / 4 + 1; i++)
            {
                i1 = i + i - 2;
                i2 = i1 + 1;
                i3 = tnn + 1 - i2;
                i4 = i3 + 1;
                wrs = twr;
                wis = twi;
                h1r = c1 * (a[i1] + a[i3]);
                h1i = c1 * (a[i2] - a[i4]);
                h2r = -(c2 * (a[i2] + a[i4]));
                h2i = c2 * (a[i1] - a[i3]);
                a[i1] = h1r + wrs * h2r - wis * h2i;
                a[i2] = h1i + wrs * h2i + wis * h2r;
                a[i3] = h1r - wrs * h2r + wis * h2i;
                a[i4] = -h1i + wrs * h2i + wis * h2r;
                twtemp = twr;
                twr = twr * twpr - twi * twpi + twr;
                twi = twi * twpr + twtemp * twpi + twi;
            }
            h1r = a[0];
            a[0] = h1r + a[1];
            a[1] = h1r - a[1];
        }

        #region ByDefinition
        /*   public static void realfftbydef(ref double[] a,
        int nn,
        bool inversefft)
    {
        double[] tmp = new double[0];
        int i = 0;
        int j = 0;
        int k = 0;
        double hre = 0;
        double him = 0;
        double re = 0;
        double im = 0;

        tmp = new double[2 * nn - 1 + 1];
        if( !inversefft )
        {
            for (i = 0; i <= nn - 1; i++)
            {
                hre = 0;
                him = 0;
                for (k = 0; k <= nn - 1; k++)
                {
                    re = a[k];
                    hre = hre + Math.Cos(2 * Math.PI * k * i / nn) * re;
                    him = him + Math.Sin(2 * Math.PI * k * i / nn) * re;
                }
                tmp[2 * i] = hre;
                tmp[2 * i + 1] = him;
            }
            for (i = 2; i <= nn - 1; i++)
            {
                a[i] = tmp[i];
            }
            a[0] = tmp[0];
            a[1] = tmp[nn];
        }
        else
        {
            tmp[0] = a[0];
            tmp[1] = 0;
            for (i = 1; i <= nn / 2 - 1; i++)
            {
                tmp[2 * i] = a[2 * i];
                tmp[2 * i + 1] = a[2 * i + 1];
                tmp[2 * nn - 2 * i] = a[2 * i];
                tmp[2 * nn - 2 * i + 1] = -a[2 * i + 1];
            }
            tmp[nn] = a[1];
            tmp[nn+1] = 0;
            for (k = 0; k <= nn - 1; k++)
            {
                hre = 0;
                him = 0;
                for (i = 0; i <= nn - 1; i++)
                {
                    re = tmp[2 * i];
                    im = tmp[2 * i + 1];
                    hre = hre + Math.Cos(-(2 * Math.PI * k * i / nn)) * re - Math.Sin(-(2 * Math.PI * k * i / nn)) * im;
                    him = him + Math.Cos(-(2 * Math.PI * k * i / nn)) * im + Math.Sin(-(2 * Math.PI * k * i / nn)) * re;
                }
                a[k] = hre / nn;
            }
        }
    }
  */
        #endregion
    }

    /// <summary>
    /// Класс реализует математические вычисления с высокой точностью.
    /// Это необходимо, например, при многоитерационном прибавлении вещественного шага к md.
    /// Вскоре начинают появляться хвосты из девяток и восьмерок.
    /// Чтобы этого не было реализован настоящий класс.
    /// </summary>
    public static class FloatHightPrecision
    {
        /// <summary>
        /// Сложение двух вещественных чисел.
        /// </summary>
        /// <param name="value1"></param>
        /// <param name="value2"></param>
        /// <returns></returns>
        public static Single Add(Single value1, Single value2)
        {
            // пока реализация самая примитивная, потом переделаю (а то как всегда все срочно надо)
            /*String vs1 = value1.ToString();
            String vs2 = value2.ToString();
            String valstr = (Single.Parse(vs1) + Single.Parse(vs2)).ToString();
            return Single.Parse(valstr);*/
            String delimiter = Thread.CurrentThread.CurrentCulture.NumberFormat.NumberDecimalSeparator;

            //String s1 = value1.ToString("F6");
            //String s1 = String.Format("{0:0.######}", Math.Round(value1, 6, MidpointRounding.ToEven));
            String s1 = Math.Round(value1, 6, MidpointRounding.ToEven).ToString("F6");
            int beg1 = 0; // Позиция начала (0 - если нет знака, 1 - если у числа есть знак)
            int sign1 = 1; // Знак числа
            if (Char.Equals(s1[0], '-'))
            {
                beg1 = 1;
                sign1 = -1;
            }
            int pointPos1 = s1.Contains(delimiter) ? s1.IndexOf(delimiter) : s1.Length;

            //String s2 = value2.ToString("F6");
            //String s2 = String.Format("{0:0.######}", Math.Round(value2, 6, MidpointRounding.ToEven));
            String s2 = Math.Round(value2, 6, MidpointRounding.ToEven).ToString("F6");
            int beg2 = 0;
            int sign2 = 1; // Знак числа
            if (Char.Equals(s2[0], '-'))
            {
                beg2 = 1;
                sign2 = -1;
            }
            int pointPos2 = s2.Contains(delimiter) ? s2.IndexOf(delimiter) : s2.Length;

            // Максимальная длина целой части
            int maxWhole = (pointPos1 - beg1) > (pointPos2 - beg2) ? (pointPos1 - beg1) : (pointPos2 - beg2);
            maxWhole++; // Для одного старшего разряда (если произойдет переход в старший разряд)
            // Максимальная длина дробной части
            int maxFraction = (s1.Length - pointPos1 - 1) > (s2.Length - pointPos2 - 1) ? (s1.Length - pointPos1 - 1) : (s2.Length - pointPos2 - 1);
            if (maxFraction < 0) maxFraction = 0;

            int[] number1 = new int[1 + maxWhole + maxFraction];
            for (int i = maxWhole, j = pointPos1 - 1; i >= 0; i--, j--)
            {
                if (j >= beg1) number1[i] = sign1 * int.Parse(s1[j].ToString());
                else number1[i] = 0;
            }
            for (int i = maxWhole + 1, j = pointPos1 + 1; i < number1.Length; i++, j++)
            {
                if (j < s1.Length) number1[i] = sign1 * int.Parse(s1[j].ToString());
                else number1[i] = 0;
            }

            int[] number2 = new int[1 + maxWhole + maxFraction];
            for (int i = maxWhole, j = pointPos2 - 1; i >= 0; i--, j--)
            {
                if (j >= beg2) number2[i] = sign2 * int.Parse(s2[j].ToString());
                else number2[i] = 0;
            }
            for (int i = maxWhole + 1, j = pointPos2 + 1; i < number2.Length; i++, j++)
            {
                if (j < s2.Length) number2[i] = sign2 * int.Parse(s2[j].ToString());
                else number2[i] = 0;
            }

            int[] result = new int[1 + maxWhole + maxFraction];
            int sign = 0;
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = number1[i] + number2[i];
                if (sign == 0 && result[i] > 0) sign = 1;
                if (sign == 0 && result[i] < 0) sign = -1;
            }

            bool going;
            do
            {
                going = false;
                for (int i = 1; i < result.Length; i++)
                {
                    if (result[i] > 9)
                    {
                        result[i - 1] += 1;
                        result[i] -= 10;
                        going = true;
                    }
                    else if (result[i] < -9)
                    {
                        result[i - 1] += -1;
                        result[i] -= 10;
                        going = true;
                    }
                    else
                    {
                        if (sign > 0 && result[i] < 0)
                        {
                            result[i - 1] += -1;
                            result[i] += 10;
                            going = true;
                        }
                        else if (sign < 0 && result[i] > 0)
                        {
                            result[i - 1] += 1;
                            result[i] += -10;
                            going = true;
                        }
                    }
                }
            } while (going);
            StringBuilder sbRes = new StringBuilder(result.Length);

            for (int i = 0; i < result.Length; i++)
                sbRes.Append(result[i] * sign);
            if (maxFraction > 0)
                sbRes.Insert(maxWhole + 1, delimiter);
            while (sbRes.Length > 0 && sbRes[0] == '0') sbRes.Remove(0, 1);
            while (sbRes.Length > 0 && sbRes[sbRes.Length - 1] == '0') sbRes.Remove(sbRes.Length - 1, 1);
            if (sbRes.Length == 1 && sbRes.ToString().Equals(delimiter))
                sbRes[0] = '0';
            if (sbRes.Length == 0) sbRes.Append('0');
            if (sign < 0) sbRes.Insert(0, '-');

            return Single.Parse(sbRes.ToString());
        }
    }
}
