using System;
using System.Collections.Generic;
using TimeZYX.Model.Spatial;
using TimeZYX.Core.Base;
using System.Collections.ObjectModel;
using TimeZYX.Model.Utils;

namespace TimeZYX.Model.Geology.Pillars
{
    // ТЗ_1_7_3_создание_модели_разломов, стр. 7
    public enum PillarKind
    {
        /// <summary>
        /// TP=0 – вертикальный двуточечный
        /// </summary>
        VERTICAL,

        /// <summary>
        /// TP=1 – наклонный двуточечный
        /// </summary>
        LINEAR,

        /// <summary>
        /// TP=2 – трехточечный
        /// </summary>
        LISTRIC,

        /// <summary>
        /// TP=3 - пятиточечный
        /// </summary>
        CURVED
    }



    //[Serializable]
    ////[DescriptionLocalized("Pilar", typeof(TimeZYX.Model.Properties.Resources))]
    //public abstract class BasePillar : DataList<PrecisedVertex>, IMoveable, IMinMax<PrecisedVertex>
    //{
        
    //    //public BasePillar(Project project, string title)
    //    //    : base(project, title)
    //    //{

    //    //}
        
    //    //public BasePillar(Project project, string title, int capacity)
    //    //    : base(project, title, capacity)
    //    //{

    //    //}

    //    public BasePillar(Project project, string title, IEnumerable<PrecisedVertex> collection, PillarKind pillarKind)
    //        : base(project, title, collection)
    //    {
    //        this.kind = pillarKind;

    //        //if (!CheckPillarKind())
    //        //    Debug.Assert(false, "Ошибка: НЕкоректный массив точек для данного типа пиллара!");

    //        //ConvertPillar(pillarKind);
    //        RecalcMinMax();
    //    }

    //    public static int GetPointCountInPillarType(PillarKind kind)
    //    {
    //        switch (kind)
    //        {
    //            case PillarKind.VERTICAL:
    //            case PillarKind.LINEAR:
    //                return 2;
    //            case PillarKind.LISTRIC:
    //                return 3;
    //            case PillarKind.CURVED:
    //                return 5;
    //            default:
    //                return 5;
    //        }
    //    }
    //    public PillarKind Kind
    //    {
    //        get
    //        {
    //            return kind;
    //        }

    //        set
    //        {
    //            ConvertPillar(value);
    //        }
    //    }

    //    public bool Attached
    //    {
    //        get
    //        {
    //            return attached;
    //        }
    //        set
    //        {
    //            attached = value;
    //        }
    //    }

    //    #region Перегруженные методы изменения списка точек в пиларе

    //    public override void Add(PrecisedVertex item)
    //    {
    //        throw new Exception("Этот метод не определен для пиларов.");

    //        //base.Add(item);
    //        //RecalcMinMax();
    //    }

    //    public override void AddRange(IEnumerable<PrecisedVertex> collection)
    //    {
    //        base.AddRange(collection);
    //        RecalcMinMax();
    //        //return;
    //        //throw new Exception("Этот метод не определен для пиларов.");
    //    }

    //    public override bool Remove(PrecisedVertex item)
    //    {
    //        //bool res = base.Remove(item);
    //        //RecalcMinMax();
    //        //return res;
    //        throw new Exception("Этот метод не определен для пиларов.");


    //    }

    //    public override void RemoveAt(int index)
    //    {
    //        //base.RemoveAt(index);
    //        //RecalcMinMax();
    //        throw new Exception("Этот метод не определен для пиларов.");
    //    }

    //    public override void Insert(int index, PrecisedVertex item)
    //    {
    //        //base.Insert(index, item);
    //        //RecalcMinMax();
    //        throw new Exception("Этот метод не определен для пиларов.");
    //    }

    //    #endregion

    //    private void RecalcMinMax()
    //    {
    //        double xmin = double.MaxValue, ymin = double.MaxValue, zmin = double.MaxValue;
    //        double xmax = double.MinValue, ymax = double.MinValue, zmax = double.MinValue;

    //        foreach (PrecisedVertex item in items)
    //        {
    //            if (xmin > item.X)
    //                xmin = item.X;
    //            if (ymin > item.Y)
    //                ymin = item.Y;
    //            if (zmin > item.Z)
    //                zmin = item.Z;

    //            if (xmax < item.X)
    //                xmax = item.X;
    //            if (ymax < item.Y)
    //                ymax = item.Y;
    //            if (zmax < item.Z)
    //                zmax = item.Z;

    //        }

    //        currentMin = new PrecisedVertex(xmin, ymin, zmin);
    //        currentMax = new PrecisedVertex(xmax, ymax, zmax);

    //    }

    //    protected virtual bool CheckPillarKind()
    //    {
    //        if (items.Count == 5)
    //            return true;

    //        return false;
    //    }

    //    protected virtual void SmoothLine()
    //    { 
        
    //    }


    //    public virtual PrecisedVertex GetCenterVertex()
    //    {
    //        PrecisedVertex v1 = items[0];
    //        PrecisedVertex v2 = items[items.Count-1];

    //        double xn = (v1.X + v2.X) / 2;
    //        double yn = (v1.Y + v2.Y) / 2;
    //        double zn = (v1.Z + v2.Z) / 2;

    //        return new PrecisedVertex(xn, yn, zn);
    //    }

    //    public virtual ReadOnlyCollection<PrecisedVertex> GetActiveVertexes()
    //    {
    //        return new ReadOnlyCollection<PrecisedVertex>(items);
    //    }

    //    #region Вспомогательные фукции пересчета точек 
    //    public virtual PrecisedVertex GetCenterVertex(PrecisedVertex v1, PrecisedVertex v2)
    //    {
    //        double xn = (v1.X + v2.X) / 2;
    //        double yn = (v1.Y + v2.Y) / 2;
    //        double zn = (v1.Z + v2.Z) / 2;

    //        return new PrecisedVertex(xn, yn, zn);
    //    }
    //    protected static void SmoothLineBy2Dots(List<PrecisedVertex> items)
    //    {
    //        if (items.Count == 5)
    //        {

    //            PrecisedVertex v1 = items[0];
    //            PrecisedVertex v2 = items[2];
    //            PrecisedVertex v3 = items[4];


    //            // расчет второй точки
    //            double x2 = (v1.X + v2.X) / 2;
    //            double y2 = (v1.Y + v2.Y) / 2;
    //            double z2 = (v1.Z + v2.Z) / 2;

    //            items[1].Move(x2, y2, z2);

    //            // расчет четвертой точки
    //            double x4 = (v2.X + v3.X) / 2;
    //            double y4 = (v2.Y + v3.Y) / 2;
    //            double z4 = (v2.Z + v3.Z) / 2;

    //            items[1].Move(x4, y4, z4);
    //        }

    //    }
    //    protected static void VerticalSmoothLine(List<PrecisedVertex> items)
    //    {
    //        if (items.Count == 5)
    //        {

    //            PrecisedVertex v1 = items[0];
    //            PrecisedVertex v2 = items[4];


    //            // расчет центральной точки
    //            double xnC = (v1.X + v2.X) / 2;
    //            double ynC = (v1.Y + v2.Y) / 2;
    //            double znC = (v1.Z + v2.Z) / 2;

    //            items[2].Move(xnC, ynC, znC);

    //            // расчет второй точки
    //            double x2 = (v1.X + xnC) / 2;
    //            double y2 = (v1.Y + ynC) / 2;
    //            double z2 = (v1.Z + znC) / 2;

    //            items[1].Move(x2, y2, z2);

    //            // расчет четвертой точки
    //            double x4 = (v2.X + xnC) / 2;
    //            double y4 = (v2.Y + ynC) / 2;
    //            double z4 = (v2.Z + znC) / 2;

    //            items[1].Move(x4, y4, z4);
    //        }

    //    }
    //    #endregion

    //    #region IMoveable Members

    //    public void Move(double dx, double dy, double dz)
    //    {
    //        foreach (var item in items)
    //        {
    //            item.Move(dx, dy, dz);
    //        }

    //        RecalcMinMax();
    //    }

    //    #endregion

    //    #region IMinMax<PrecisedVertex> Members

    //    public PrecisedVertex Min
    //    {
    //        get 
    //        { 
    //            return currentMin;
    //        }
    //    }

    //    public PrecisedVertex Max
    //    {
    //        get 
    //        { 
    //            return currentMax;
    //        }
    //    }

    //    #endregion

    //    #region Вспомогательные функции

    //    /// <summary>
    //    /// Преобразование пилара из одного типа в другой, или преобразовать массив точек (items) 
    //    /// к указанному типу пилара.
    //    /// </summary>
    //    /// <param name="newKind">Тип пилара</param>
    //    private void ConvertPillar(PillarKind newKind)
    //    {
    //        PrecisedVertex[] buf=null;
    //        PrecisedVertex center;
    //        switch (newKind)
    //        {
    //            case PillarKind.VERTICAL:
    //                buf=new PrecisedVertex[2];
    //                center = GetCenterVertex();
    //                buf[0]=new PrecisedVertex(center.X,center.Y,items[0].Z);
    //                buf[1] = new PrecisedVertex(center.X, center.Y, items[items.Count - 1].Z);
    //                break;
    //            case PillarKind.LINEAR:
    //                buf = new PrecisedVertex[2];
    //                buf[0] = items[0];
    //                buf[1] = items[items.Count - 1];
    //                break;
    //            case PillarKind.CURVED:
    //                buf = new PrecisedVertex[3];
    //                if (items.Count < 3)
    //                {
    //                    buf[0] = items[0];
    //                    buf[2] = items[items.Count - 1];
    //                    buf[1] = GetCenterVertex(buf[0], buf[2]);
    //                    Clear();
    //                    AddRange(buf);
    //                }
    //                else
    //                {
    //                    buf[0] = items[0];
    //                    buf[1] = GetCenterVertex();
    //                    buf[2] = items[items.Count - 1];
    //                }
    //                break;
    //            case PillarKind.LISTRIC:
    //                buf = new PrecisedVertex[5];
    //                if (items.Count < 2)
    //                {
    //                    buf[0] = items[0];
    //                    buf[2] = GetCenterVertex(items[0], items[1]);
    //                    buf[1] = GetCenterVertex(buf[0], buf[1]);
    //                    buf[3] = GetCenterVertex(buf[2], buf[4]);
    //                    buf[4] = items[items.Count - 1];
    //                }
    //                else
    //                    if (items.Count == 3)
    //                    {
    //                        buf[0] = items[0];
    //                        buf[2] = items[1];
    //                        buf[1] = GetCenterVertex(buf[0], buf[1]);
    //                        buf[3] = GetCenterVertex(buf[2], buf[4]);
    //                        buf[4] = items[2];
    //                    }
    //                    else 
    //                    buf = items.ToArray();
    //                break;
    //        }
    //        items =new List<PrecisedVertex>( buf);
    //        //Clear();
    //        //AddRange(items);
    //        return;

    //    }

    //    public bool SetPillarPoint(int numPoint, PrecisedVertex newLoaction)
    //    {
    //        try
    //        {
    //            if (numPoint < 0 || numPoint > Count - 1)
    //                return false;
    //            int minIndex = -1;
    //            int maxIndex = -1;
    //            if (numPoint == 0)
    //            {
    //                maxIndex = 1;
    //            }
    //            else if (numPoint == Count - 1)
    //            {
    //                minIndex = Count - 2;
    //            }
    //            else
    //            {
    //                maxIndex = numPoint + 1;
    //                minIndex = numPoint - 1;
    //            }

    //            if (minIndex != -1 && maxIndex != -1)
    //            {
    //                if (newLoaction.Z >= items[maxIndex].Z
    //                   || newLoaction.Z <= items[minIndex].Z)
    //                    return false;
    //                else
    //                    items[numPoint] = newLoaction;
    //            }
    //            else if (minIndex != -1)
    //            {
    //                if (newLoaction.Z <= items[minIndex].Z)
    //                    return false;
    //                else
    //                    items[numPoint] = newLoaction;
    //            }
    //            else
    //            {
    //                if (newLoaction.Z >= items[maxIndex].Z)
    //                    return false;
    //                else
    //                    items[numPoint] = newLoaction;
    //            }

    //            return true;
    //        }
    //        catch
    //        {
    //            return false;
    //        }
    //    }


    //    /// <summary>
    //    /// возвращает точку на параметризованной кривой
    //    /// параметр от 0 до 1. Кривая параметризируется в 
    //    /// соответствии с координатой z.
    //    /// 
    //    /// Нужно уточнить нужна ли эта функция вне этого класса?
    //    /// </summary>
    //    /// <param name="t"></param>
    //    /// <returns></returns>
    //    public PrecisedVertex ParametricByZ(double z, PrecisedVertex[] attractors)
    //    {

    //        // заглушка
    //        return ParametricByT(CrossPoint(z, attractors));
    //    }
    //    public PrecisedVertex ParametricByZ(double z)
    //    {

    //        // заглушка
    //        return ParametricByT(CrossPoint(z));
    //    }

    //    /// <summary>
    //    /// возвращает точку на параметризованной кривой
    //    /// параметр от 0 до 1. Кривая параметризируется в 
    //    /// соответствии с координатой z
    //    /// 
    //    /// Нужно уточнить нужна ли эта функция вне этого класса?
    //    /// </summary>
    //    /// <param name="t"></param>
    //    /// <returns></returns>
    //    public virtual PrecisedVertex ParametricByT(double t, PrecisedVertex[] attractors)
    //    {
            
    //        return new PrecisedVertex();
    //    }

    //    public virtual PrecisedVertex ParametricByT(double t)
    //    {
            
    //        return new PrecisedVertex();
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
    //    public virtual double CrossPoint(double z)
    //    {
    //        return 0;
    //    }

    //    public virtual double CrossPoint(double z,PrecisedVertex[] attractors)
    //    {
    //        return 0;
    //    }

    //    /// <summary>
    //    /// Перемещение пиллара с его притяжками
    //    /// </summary>
    //    /// <param name="x"></param>
    //    /// <param name="y"></param>
    //    /// <param name="z"></param>
    //    /// <returns></returns>
    //    public  void Move(double dx, double dy, double dz, PrecisedVertex[] attractors)
    //    {
    //        Move(dx, dy, dz);

    //        foreach (var attractor in attractors)
    //        {
    //            attractor.Move(dx, dy, dz);
    //        }

    //        return;
    //    }
    //    /// <summary>
    //    /// Возможно тут вместо  double dx, double dy, double x, double y 
    //    /// надо передавать PrecisedVertex
    //    /// </summary>
    //    /// <param name="surface"></param>
    //    /// <param name="dx"></param>
    //    /// <param name="dy"></param>
    //    /// <param name="x"></param>
    //    /// <param name="y"></param>
    //    /// <returns></returns>
    //     double PointOnSurface(Surface surface, double x, double y)
    //    {

    //        double z = -999;

    //        int i0 = (int)Math.Floor((x - surface.Xmin) / surface.XStep);
    //        int i1 = i0 + 1;
    //        // int ix = (Math.Abs(x - xmin2 - dx2 * c) < Math.Abs(x - xmin2 - dx2 * f)) ? c : f;

    //        if (i0 < 0 || i1 > (surface.Nx - 1)) { return z; }



    //        int j0 = (int)Math.Floor((y - surface.Ymin) / surface.YStep);
    //        int j1 = j0 + 1;
    //        // int iy= (Math.Abs(y - ymin2 - dy2 * c) < Math.Abs(y - ymin2 - dy2 * f)) ? c : f;

    //        if (j0 < 0 || j1 > (surface.Ny - 1)) { return z; }

    //        if (surface[i0, j0].Z == surface.NullValue || surface[i1, j0].Z == surface.NullValue || surface[i1, j1].Z == surface.NullValue || surface[i0, j1].Z == surface.NullValue)
    //            return z;
    //        double xx = (x - surface.Xmin - i0 * surface.XStep) / surface.XStep;
    //        double yy = (y - surface.Ymin - j0 * surface.YStep) / surface.YStep;
    //        double nx, ny, nz;
    //        if (yy - xx > 0)
    //        {

    //            nx = (j0 - j1) * (surface[i1, j1].Z - surface[i0, j1].Z);
    //            ny = (surface[i0, j0].Z - surface[i0, j1].Z) * (i1 - i0);
    //            nz = -(j0 - j1) * (i1 - i0);

    //            z = (nx * (0 - xx) + ny * (1 - yy)) / nz + surface[i0, j1].Z;
    //        }
    //        else
    //        {
    //            nx = -(j1 - j0) * (surface[i0, j0].Z - surface[i1, j0].Z);
    //            ny = -(surface[i1, j1].Z - surface[i1, j0].Z) * (i0 - i1);
    //            nz = (j1 - j0) * (i0 - i1);

    //            z = (nx * (1 - xx) + ny * (0 - yy)) / nz + surface[i1, j0].Z;
    //        }

    //        return z;
    //    }

    //    /// <summary>
    //    /// Возможно тут вместо  double dx, double dy, double x, double y 
    //    /// надо передавать PrecisedVertex
    //    /// Скоее всего эта функция должна быть private ???
    //    /// 
    //    /// </summary>
    //    /// <param name="surface"></param>
    //    /// <param name="dx"></param>
    //    /// <param name="dy"></param>
    //    /// <param name="x"></param>
    //    /// <param name="y"></param>
    //    /// <returns></returns>
    //     double[] CrossSur(Surface surface,  double t1, double t2,out double t0)
    //    {
    //        t0 = -999;
    //        double[] res = new double[3];
    //        double[] p1 = new double[3];
    //        double[] p2 = new double[3];
    //        PrecisedVertex pnt = this.ParametricByT(t1);
    //        p1[0] = pnt.X;
    //        p1[1] = pnt.Y;
    //        p1[2] = pnt.Z;
    //        pnt = this.ParametricByT(t2);
    //        p2[0] = pnt.X;
    //        p2[1] = pnt.Y;
    //        p2[2] = pnt.Z;
    //        double z1 = PointOnSurface(surface, p1[0], p1[1]);
    //        double z2 = PointOnSurface(surface, p2[0], p2[1]);
    //        if (z1 == -999 || z2 == -999) return null;
    //        if ((p1[2] - z1) * (p1[2] - z1) < 0.0001)
    //        {
    //            t0 = t1;
    //            return p1;
    //        }
    //        else if ((p2[2] - z2) * (p2[2] - z2) < 0.0001)
    //        {
    //            t0 = t2;
    //            return p2;
    //        }

    //        else if ((p1[2] - p2[2]) * (p1[2] - p2[2]) < 0.0001)
    //        {
    //            res[0] = (p1[0] + p2[0]) * 0.5;
    //            res[1] = (p1[1] + p2[1]) * 0.5;
    //            res[2] = PointOnSurface(surface, res[0], res[1]);
    //            if (res[2] == -999) return null;
    //            t0 = (t1 + t2) * 0.5;
    //            return res;
    //        }
    //        double[] p3 = new double[3];
    //        pnt = this.ParametricByT((t1 + t2) * 0.5);
    //        p3[0] = pnt.X;
    //        p3[1] = pnt.Y;
    //        p3[2] = pnt.Z;
    //        double z3 = PointOnSurface(surface, p3[0], p3[1]);
    //        if (z3 == -999) return null;
    //        if ((p3[2] - z3) * (p3[2] - z3) < 0.0001)
    //        {
    //            t0 = (t1 + t2) * 0.5;
    //            return p3;
    //        }
    //        else if (p3[2] > z3)
    //            return CrossSur(surface, t1, (t2 + t1) * 0.5, out t0);
    //        else if (p3[2] < z3)
    //            return CrossSur(surface, (t2 + t1) * 0.5, t2, out t0);

    //        else
    //        {

    //            return null;
    //        }
    //    }
    //     double[] CrossSur(Surface surface, double t1, double t2,out double t0, PrecisedVertex[] attractors)
    //     {
    //         t0 = -999;
    //         double[] res = new double[3];
    //         double[] p1 = new double[3];
    //         double[] p2 = new double[3];
    //         PrecisedVertex pnt = this.ParametricByT(t1, attractors);
    //         p1[0] = pnt.X;
    //         p1[1] = pnt.Y;
    //         p1[2] = pnt.Z;
    //         pnt = this.ParametricByT(t2, attractors);
    //         p2[0] = pnt.X;
    //         p2[1] = pnt.Y;
    //         p2[2] = pnt.Z;
    //         double z1 = PointOnSurface(surface, p1[0], p1[1]);
    //         double z2 = PointOnSurface(surface, p2[0], p2[1]);
    //         if (z1 == -999 || z2 == -999) return null;
    //         if ((p1[2] - z1) * (p1[2] - z1) < 0.0001)
    //         {
    //             t0 = t1;
    //             return p1;
    //         }
    //         else if ((p2[2] - z2) * (p2[2] - z2) < 0.0001)
    //         {
    //             t0 = t2;
    //             return p2;
    //         }

    //         else if ((p1[2] - p2[2]) * (p1[2] - p2[2]) < 0.0001)
    //         {
    //             res[0] = (p1[0] + p2[0]) * 0.5;
    //             res[1] = (p1[1] + p2[1]) * 0.5;
    //             res[2] = PointOnSurface(surface, res[0], res[1]);
    //             if (res[2] == -999) return null;
    //             t0 = (t1 + t2) * 0.5;
    //             return res;
    //         }
    //         double[] p3 = new double[3];
    //         pnt = this.ParametricByT((t1 + t2) * 0.5, attractors);
    //         p3[0] = pnt.X;
    //         p3[1] = pnt.Y;
    //         p3[2] = pnt.Z;
    //         double z3 = PointOnSurface(surface, p3[0], p3[1]);
    //         if (z3 == -999) return null;
    //         if ((p3[2] - z3) * (p3[2] - z3) < 0.0001)
    //         {
    //             t0 = (t1 + t2) * 0.5;
    //             return p3;
    //         }
    //         else if (p3[2] > z3)
    //             return CrossSur(surface, t1, (t2 + t1) * 0.5, out t0,attractors);
    //         else if (p3[2] < z3)
    //             return CrossSur(surface, (t2 + t1) * 0.5, t2, out t0, attractors);

    //         else
    //         {

    //             return null;
    //         }
    //     }

    //    /// <summary>
    //    /// 
    //    /// </summary>
    //    /// <param name="surface"></param>
    //    /// <param name="dx"></param>
    //    /// <param name="dy"></param>
    //    /// <returns></returns>
    //    public void CrossSurface(Surface surface, out double z,out double t)
    //    {
    //        t = -999;
    //        z = -999;
    //        double t0 = CrossPoint(surface.Min.Z);
    //        double t1 = CrossPoint(surface.Max.Z);
    //        double[] res = CrossSur(surface, t0, t1, out t);
    //        if (res != null)
    //        {

    //            z = ParametricByT(t).Z;
    //        }
    //        return;
    //    }
    //    public void CrossSurface(Surface surface, out double z, out double t,PrecisedVertex[] attractors)
    //    {
    //        t = -999;
    //        z = -999;
    //        double t0 = CrossPoint(surface.Min.Z);
    //        double t1 = CrossPoint(surface.Max.Z);
    //        double[] res = CrossSur(surface, t0, t1, out t, attractors);
    //        if (res != null)
    //        {

    //            z = ParametricByT(t, attractors).Z;
    //        }
    //        return;
    //    }
       
    //    #endregion

    //    #region Создание точек притяжки для отрисовки и алгоритмов
    //    public virtual PrecisedVertex[] GetBesierPoints()
    //    {
    //        int BezierParam = 10;
    //        PrecisedVertex[] ctrlpointsL;
    //        PrecisedVertex[] ctrlpointsR;
    //        MathFunctions.GetControlPoints(items.ToArray(), out ctrlpointsL, out ctrlpointsR);

    //        double t;
    //        int iter = 0;
    //        PrecisedVertex[] bezierPoints = new PrecisedVertex[(Count - 1) * (BezierParam+1)];
    //        for (int i = 0; i < items.Count - 1; i++)
    //        {
    //            for (int j = 0; j <= BezierParam; j++)
    //            {
    //                t = j * 1.0 / BezierParam;
    //                PrecisedVertex V = MathFunctions.BezierCurve(items[i], ctrlpointsL[i], ctrlpointsR[i], items[i + 1], t);
    //                bezierPoints[iter++] = V;
    //            }
    //        }
    //        return null;
    //    }
    //    #endregion
    //    PrecisedVertex currentMin;
    //    PrecisedVertex currentMax;
    //    PillarKind kind = PillarKind.VERTICAL;
    //    bool attached = false;
    //}
}
