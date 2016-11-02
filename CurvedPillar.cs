using System;
using System.Collections.Generic;
using System.Text;
using TimeZYX.Core.Base;
using TimeZYX.Model.Spatial;

namespace TimeZYX.Model.Geology.Pillars
{
    //[Serializable]
    //public class CurvedPillar : BasePillar
    //{
    //    public CurvedPillar(Project project, string title, IEnumerable<PrecisedVertex> collection)
    //        : base(project, title, collection, PillarKind.CURVED)
    //    { 
        
    //    }

    //    protected override bool CheckPillarKind()
    //    {
    //        if (items.Count == 2)
    //        {
    //            items.Insert(1, new PrecisedVertex());
    //            items.Insert(1, new PrecisedVertex());
    //            items.Insert(1, new PrecisedVertex());

    //            VerticalSmoothLine(items);
    //            return true;
    //        }
    //        else if (items.Count == 3)
    //        {
    //            items.Insert(1, new PrecisedVertex());
    //            items.Insert(3, new PrecisedVertex());

    //            SmoothLineBy2Dots(items);

    //            return true;

    //        }
    //        else if (items.Count == 5)
    //        {
    //            return true;
    //        }
            
    //        return false;
    //    }

    //    #region Создание точек притяжки для отрисовки и алгоритмов
    //    public override PrecisedVertex[] GetBesierPoints()
    //    {
    //        PrecisedVertex[] attractors = new PrecisedVertex[items.Count];
    //        PrecisedVertex tmp, att;
    //        tmp = new PrecisedVertex(items[2].X - items[0].X, items[2].Y - items[0].Y, items[2].Z - items[0].Z);
    //        double l = Math.Sqrt(tmp.X * tmp.X + tmp.Y * tmp.Y + tmp.Z * tmp.Z);
    //        double l1 = (tmp.X * (items[1].X - items[0].X) +
    //        tmp.Y * (items[1].Y - items[0].Y) +
    //        tmp.Z * (items[1].Z - items[0].Z)) / l;

    //        double l2 = (tmp.X * (items[2].X - items[1].X) +
    //        tmp.Y * (items[2].Y - items[1].Y) +
    //        tmp.Z * (items[2].Z - items[1].Z)) / l;
    //        tmp = new PrecisedVertex(tmp.X / l / 5, tmp.Y / l / 5, tmp.Z / l / 5);

    //        attractors[0] = new PrecisedVertex(items[1].X - tmp.X * l1, items[1].Y - tmp.Y * l1, items[1].Z - tmp.Z * l1);
    //        attractors[1] = new PrecisedVertex(items[1].X + tmp.X * l2, items[1].Y + tmp.Y * l2, items[1].Z + tmp.Z * l2);
    //        //заглушка
    //        return attractors;
    //    }

    //    #endregion
    //     #region внутренние функции
    //    double[] QuadricRoot(double b, double c, double d)
    //    {
    //        double[] troot;
    //        if (b == 0)
    //        {
    //            if (c == 0)
    //            {
    //                troot = new double[0];
    //                return troot;
    //            }
    //            troot = new double[] { -d / c };
    //            return troot;
    //        }
    //       double  des = c * c - b * d * 4;
    //        if (des < 0)
    //        {
    //            troot = new double[0];
    //            return troot;
    //        }
    //        if (des == 0)
    //        {
    //            troot = new double[] { -c / (2 * b) };
    //            return troot;
    //        }

    //        troot = new double[] { (-c - Math.Sqrt(des)) / (2 * b), (-c + Math.Sqrt(des)) / (2 * b) };
    //        return troot;

    //    }
    //    #endregion
    //    #region Вспомогательные функции
    //    /// <summary>
    //    /// возвращает точку на параметризованной кривой
    //    /// параметр от 0 до 1. Кривая параметризируется в 
    //    /// соответствии с координатой z
    //    /// 
    //    /// Нужно уточнить нужна ли эта функция вне этого класса?
    //    /// </summary>
    //    /// <param name="t"></param>
    //    /// <returns></returns>
    //    public override PrecisedVertex ParametricByT(double t, PrecisedVertex[] attractors)
    //    {
    //        double[] res = new double[3]; 
    //        if (t <= 0)
    //        {
    //            res[0] = items[0].X + 2 * (attractors[0].X - items[0].X) * t;
    //            res[1] = items[0].Y + 2 * (attractors[0].Y - items[0].Y) * t;
    //            res[2] = items[0].Z + 2 * (attractors[0].Z - items[0].Z) * t;
    //        }
    //        else if (t >= 1)
    //        {
    //            res[0] = items[2].X - 2 * (attractors[1].X - items[2].X) * (t - 1);
    //            res[1] = items[2].Y - 2 * (attractors[1].Y - items[2].Y) * (t - 1);
    //            res[2] = items[2].Z - 2 * (attractors[1].Z - items[2].Z) * (t - 1);
    //        }
    //        else if (t < 0.5)
    //        {
    //            t = t * 2;
    //            res[0] = items[0].X * (1 - t) * (1 - t) +
    //                2 * attractors[0].X * t * (1 - t) + items[1].X * t * t;
    //            res[1] = items[0].Y * (1 - t) * (1 - t) +
    //                2 * attractors[0].Y * t * (1 - t) + items[1].Y * t * t;
    //            res[2] = items[0].Z * (1 - t) * (1 - t) +
    //                2 * attractors[0].Z * t * (1 - t) + items[1].Z * t * t;
    //        }
    //        else
    //        {
    //            t = (t - 0.5) * 2;
    //            res[0] = items[1].X * (1 - t) * (1 - t) +
    //                    2 * attractors[1].X * t * (1 - t) + items[2].X * t * t;
    //            res[1] = items[1].Y * (1 - t) * (1 - t) +
    //                2 * attractors[1].Y * t * (1 - t) + items[2].Y * t * t;
    //            res[2] = items[1].Z * (1 - t) * (1 - t) +
    //                2 * attractors[1].Z * t * (1 - t) + items[2].Z * t * t;
    //        }
    //        // заглушка
    //        return new PrecisedVertex(res[0], res[1],res[2]);
    //    }

    //    public override PrecisedVertex GetCenterVertex()
    //    {
    //        return items[2];
    //    }
        
    //    /// <summary>
    //    /// возвращает параметр t, при котором пилар пересекает плоскость z
    //    /// 
    //    /// Нужно уточнить нужна ли эта функция вне этого класса?
    //    /// </summary>
    //    /// <param name="t1"></param>
    //    /// <param name="t2"></param>
    //    /// <param name="z"></param>
    //    /// <returns></returns>
    //    public override double CrossPoint(double z, PrecisedVertex[] attractors)
    //    {
    //        double[] troot;
    //        double eps=0.0001, a, b, c, d;
    //        if (z < items[0].Z + eps)
    //        {

    //            a = 2 * (attractors[0].Z - items[0].Z);
    //            b = items[0].Z - z;
    //            return -b / a;

    //        }
    //        if (z > items[items.Count - 1].Z + eps)
    //        {

    //            a = 2 * (attractors[1].Z - items[2].Z);
    //            b = items[2].Z - z;
    //            return b / a + 1;

    //        }
    //        if (z < items[1].Z + eps)
    //        {
    //            a = 0;
    //            b = items[0].Z + items[1].Z - 2 * attractors[0].Z;
    //            c = -2 * items[0].Z + 2 * attractors[0].Z;
    //            d = items[0].Z - z;
    //            troot = QuadricRoot(b, c, d);
    //            for (int t = 0; t < troot.Length; t++)
    //                if (troot[t] >= 0 - eps && troot[t] <= 1 + eps) return troot[t] * 0.5;
    //        }
    //        if (z < items[2].Z + eps)
    //        {
    //            a = 0;
    //            b = items[1].Z + items[2].Z - 2 * attractors[1].Z;
    //            c = -2 * items[1].Z + 2 * attractors[1].Z;
    //            d = items[1].Z - z;

    //            troot= QuadricRoot( b, c, d);
    //            for (int t = 0; t < troot.Length; t++)
    //                if (troot[t] >= 0 - eps && troot[t] <= 1 + eps) return troot[t] * 0.5 + 0.5;
    //        }
    //        return -999;
    //    }
    //    #endregion

    //}
}
