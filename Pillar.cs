using System;
using System.Collections.Generic;
using TimeZYX.Model.Spatial;
using TimeZYX.Core.Base;

namespace TimeZYX.Model.Geology.Pillars
{
    [Serializable]
    public class Pillar : Data, IMinMax<PrecisedVertex>, IMoveable
    {
        private readonly IList<PrecisedVertex> vertices;
        private PillarImpl impl;
        private PrecisedVertex currentMin;
        private PrecisedVertex currentMax;
       
        private PillarKind kind = PillarKind.LINEAR;
        private bool attached;

        public bool Attached
        {
            get { return attached; }
            set { attached = value; }
        }

        public Pillar(Project project, string title, IEnumerable<PrecisedVertex> items, PillarImpl impl)
            :base(project,title)
        {
            IEnumerator<PrecisedVertex> iter = items.GetEnumerator();
            int count = 0;
            while (iter.MoveNext())
                count++;
            System.Diagnostics.Debug.Assert(count >= 5);
            this.impl = impl;
           List<PrecisedVertex> verts = new List<PrecisedVertex>(5);
            foreach (var item in items)
                verts.Add(item);
            if (verts[0].Z > verts[1].Z)
                verts.Reverse();
            vertices = verts; 
               // SetVertexes(verts.ToArray());
                //impl.ConvertPillar(items);
                impl.ConvertPillar(vertices);
            RecalcMinMax();
        }
        public Pillar(Project project, string title, IEnumerable<PrecisedVertex> items)
            :this(project,title,items,new LinearPillarImpl())
        {

        }

        public PrecisedVertex Center()
        {
            return impl.GetCenter(vertices);
        }

        public PrecisedVertex GetVertex(int index)
        {
            return impl.GetVertex(vertices, index);
        }

        public void SetVertex(PrecisedVertex newVertex, int index)
        {
            impl.SetVertex(vertices, newVertex, index);
        }
       
        public void SetVertexes(PrecisedVertex[] newVertexes)
        {
            impl.SetVertexes(vertices, newVertexes);
        }

        public IList<PrecisedVertex> GetControlPoints()
        {
            return impl.GetVertices(vertices);
        }

        public IList<PrecisedVertex> GetAllPoints()
        {
            return vertices;
        }

        protected virtual void RecalcMinMax()
        {
            double xmin = double.MaxValue, ymin = double.MaxValue, zmin = double.MaxValue;
            double xmax = double.MinValue, ymax = double.MinValue, zmax = double.MinValue;

            foreach (PrecisedVertex item in vertices)
            {
                if (xmin > item.X)
                    xmin = item.X;
                if (ymin > item.Y)
                    ymin = item.Y;
                if (zmin > item.Z)
                    zmin = item.Z;

                if (xmax < item.X)
                    xmax = item.X;
                if (ymax < item.Y)
                    ymax = item.Y;
                if (zmax < item.Z)
                    zmax = item.Z;

            }

            currentMin = new PrecisedVertex(xmin, ymin, zmin);
            currentMax = new PrecisedVertex(xmax, ymax, zmax);

        }

        public PillarKind Kind
        {
            get
            {
                return kind;
            }

            set
            {
                //if (kind == value)
                //    return;
                switch (value)
                {
                    case PillarKind.VERTICAL:
                        SetImpl(new VerticalPillarImpl());
                        break;
                    case PillarKind.LISTRIC:
                        SetImpl(new ListricPillarImpl());
                        break;
                    case PillarKind.CURVED:
                        SetImpl(new CurvedPillarImpl());
                        break;
                    case PillarKind.LINEAR:
                        SetImpl(new LinearPillarImpl());
                        break;
                    default:
                        SetImpl(new LinearPillarImpl());
                        break;
                }
                kind = value;
            }
        }
       
        private void GetControlvectors(out PrecisedVertex[] P1, out PrecisedVertex[] P2, out PrecisedVertex beginVector, out PrecisedVertex endVector)
        {
            
            TimeZYX.Model.Utils.MathFunctions.GetControlPoints(new List<PrecisedVertex>(vertices).ToArray(), out P1, out P2);
            beginVector = new PrecisedVertex(-3 * vertices[0].X + 3 * P1[0].X, -3 * vertices[0].Y + 3 * P1[0].Y, (vertices[1].Z - vertices[0].Z));
            double l = Math.Sqrt(beginVector.X * beginVector.X + beginVector.Y * beginVector.Y + beginVector.Z * beginVector.Z);
            beginVector = new PrecisedVertex(beginVector.X / l, beginVector.Y / l, beginVector.Z / l);

            endVector = new PrecisedVertex(3 * vertices[4].X - 3 * P2[3].X, 3 * vertices[4].Y - 3 * P2[3].Y, (vertices[4].Z - vertices[3].Z));
            l = Math.Sqrt(endVector.X * endVector.X + endVector.Y * endVector.Y + endVector.Z * endVector.Z);
            endVector = new PrecisedVertex(endVector.X / l, endVector.Y / l, endVector.Z / l);
            
        }
       
        private double CrossSur(BaseSurface surface, double t1, double t2)
        {

            PrecisedVertex p1 = this.ParametricByT(t1);

            PrecisedVertex p2 = this.ParametricByT(t2);
            int i0 = (int)Math.Floor((p1.X - surface.Xmin) / surface.XStep);
            if (i0 >= 0 && i0 < surface.Nx - 1 && i0 == (int)Math.Floor((p2.X - surface.Xmin) / surface.XStep))
            {
                int j0 = (int)Math.Floor((p1.Y - surface.Ymin) / surface.YStep);
                if (j0 >= 0 && j0 < surface.Ny - 1 && j0 == (int)Math.Floor((p2.Y - surface.Ymin) / surface.YStep))
                {
                    int i1 = i0 + 1;
                    int j1 = j0 + 1;

                    double nx, ny, nz, ax, ay, az, t;
                    if (surface[i0, j0].Z != surface.NullValue && (surface[i1, j1].Z != surface.NullValue))
                    {
                        ax = p2.X - p1.X;
                        ay = p2.Y - p1.Y;
                        az = p2.Z - p1.Z;
                        if (surface[i0, j1].Z != surface.NullValue)
                        {
                            nx = (surface[i0, j0].Y - surface[i0, j1].Y) * (surface[i1, j1].Z - surface[i0, j1].Z);
                            ny = (surface[i0, j0].Z - surface[i0, j1].Z) * (surface[i1, j1].X - surface[i0, j1].X);
                            nz = -(surface[i0, j0].Y - surface[i0, j1].Y) * (surface[i1, j1].X - surface[i0, j1].X);

                            t = nx * (surface[i0, j0].X - p1.X) + ny * (surface[i0, j0].Y - p1.Y) + nz * (surface[i0, j0].Z - p1.Z);
                            t /= ax * nx + ay * ny + az * nz;
                            if (t >= 0 && t <= 1)
                            {
                                p2 = new PrecisedVertex(p1.X + ax * t, p1.Y + ay * t, p1.Z + az * t);
                                if ((p2.Y - surface[i0, j0].Y) / surface.YStep > (p2.X - surface[i0, j0].X) / surface.XStep)
                                    return p2.Z;
                            }
                        }
                        if (surface[i1, j0].Z != surface.NullValue)
                        {
                            nx = (surface[i1, j1].Y - surface[i1, j0].Y) * (surface[i0, j0].Z - surface[i1, j0].Z);
                            ny = (surface[i1, j1].Z - surface[i1, j0].Z) * (surface[i0, j0].X - surface[i1, j0].X);
                            nz = -(surface[i1, j1].Y - surface[i1, j0].Y) * (surface[i0, j0].X - surface[i1, j0].X);
                            t = nx * (surface[i0, j0].X - p1.X) + ny * (surface[i0, j0].Y - p1.Y) + nz * (surface[i0, j0].Z - p1.Z);
                            t /= ax * nx + ay * ny + az * nz;
                            if (t >= 0 && t <= 1)
                            {
                                p2 = new PrecisedVertex(p1.X + ax * t, p1.Y + ay * t, p1.Z + az * t);
                                if ((p2.Y - surface[i0, j0].Y) / surface.YStep <= (p2.X - surface[i0, j0].X) / surface.XStep)
                                    return p2.Z;
                            }
                        }
                    }


                }

            }

            double z1 = PointOnSurface(surface, p1.X, p1.Y);
            double z2 = PointOnSurface(surface, p2.X, p2.Y);
            double tbuf = -999;
            PrecisedVertex pbuf = new PrecisedVertex(-999, -999, -999);
            double zbuf = -999;
            if (z1 == -999 && z2 == -999)
            {

                int ind0 = (int)Math.Ceiling(t1 * 4);
                ind0 = Math.Max(0, ind0);
                int ind1 = (int)Math.Floor(t2 * 4);
                ind1 = Math.Min(4, ind1);
                for (int v = ind0; v <= ind1; v++)
                {
                    zbuf = PointOnSurface(surface, vertices[v].X, vertices[v].Y);
                    if (zbuf != -999)
                    {
                        if (zbuf < vertices[v].Z)
                        {
                            z1 = zbuf;
                            t1 = v * 0.25;
                            p1 = vertices[v];
                            break;
                        }
                        else
                        {
                            z2 = zbuf;
                            t2 = v * 0.25;
                            p2 = vertices[v];
                            break;
                        }
                    }
                }
                if (zbuf == -999) return -999;
            }
            if (z1 == -999)
            {
                do
                {
                    if (Math.Abs(t1 - t2) < 0.0001)
                    {
                        if (z1 == -999 || (p2.Z - z2) * (p1.Z - z1) > 0) return -999;
                        else break;
                    }
                    tbuf = (t1 + t2) * 0.5;
                    pbuf = this.ParametricByT(tbuf);
                    zbuf = PointOnSurface(surface, pbuf.X, pbuf.Y);
                    if (zbuf != -999)
                    {
                        if ((pbuf.Z - zbuf) * (p2.Z - z2) > 0)
                        {
                            t2 = tbuf;
                            z2 = zbuf;
                            p2 = pbuf;
                        }
                        else
                        {
                            t1 = tbuf;
                            z1 = zbuf;
                            p1 = pbuf;
                            break;
                        }
                    }
                    else
                    {
                        t1 = tbuf;
                        z1 = zbuf;
                        p1 = pbuf;
                    }
                } while ((zbuf == -999) || (zbuf - pbuf.Z) * (zbuf - pbuf.Z) > 0.0001);
                if ((zbuf - pbuf.Z) * (zbuf - pbuf.Z) <= 0.0001) return zbuf;
            }
            else if (z2 == -999)
            {
                do
                {
                    if (Math.Abs(t1 - t2) < 0.0001)
                    {
                        if (z2 == -999 || (p2.Z - z2) * (p1.Z - z1) > 0) return -999;
                        else break;
                    }
                    tbuf = (t1 + t2) * 0.5;
                    pbuf = this.ParametricByT(tbuf);
                    zbuf = PointOnSurface(surface, pbuf.X, pbuf.Y);
                    if (zbuf != -999)
                    {
                        if ((pbuf.Z - zbuf) * (p1.Z - z1) > 0)
                        {
                            t1 = tbuf;
                            z1 = zbuf;
                            p1 = pbuf;
                        }
                        else
                        {
                            t2 = tbuf;
                            z2 = zbuf;
                            p2 = pbuf;
                            break;
                        }
                    }
                    else
                    {
                        t2 = tbuf;
                        z2 = zbuf;
                        p2 = pbuf;
                    }
                } while ((zbuf == -999) || (zbuf - pbuf.Z) * (zbuf - pbuf.Z) > 0.0001);
                if ((zbuf - pbuf.Z) * (zbuf - pbuf.Z) <= 0.0001) return zbuf;
            }

            if ((p1.Z - z1) * (p1.Z - z1) < 0.0001)
            {

                return p1.Z;
            }
            else if ((p2.Z - z2) * (p2.Z - z2) < 0.0001)
            {
                return p2.Z;
            }

            else if ((p1.Z - p2.Z) * (p1.Z - p2.Z) < 0.0001)
            {
                return PointOnSurface(surface, (p1.X + p2.X) * 0.5, (p1.Y + p2.Y) * 0.5);
            }

            PrecisedVertex p3 = this.ParametricByT((t1 + t2) * 0.5);
            double z3 = PointOnSurface(surface, p3.X, p3.Y);
            if (z3 == -999) return z3;
            if ((p3.Z - z3) * (p3.Z - z3) < 0.0001)
            {
                return p3.Z;
            }
            else if (p3.Z > z3)
                return CrossSur(surface, t1, (t2 + t1) * 0.5);
            else if (p3.Z < z3)
                return CrossSur(surface, (t2 + t1) * 0.5, t2);

            else
            {

                return -999;
            }
        }

        private static PrecisedVertex BezierCurve(PrecisedVertex P0, PrecisedVertex P1, PrecisedVertex P2, PrecisedVertex P3, double t)
        {
            double x = P0.X * Math.Pow((1 - t), 3) + 3 * Math.Pow((1 - t), 2) * t * P1.X + 3 * (1 - t) * t * t * P2.X + t * t * t * P3.X;
            double y = P0.Y * Math.Pow((1 - t), 3) + 3 * Math.Pow((1 - t), 2) * t * P1.Y + 3 * (1 - t) * t * t * P2.Y + t * t * t * P3.Y;
            double z = P0.Z + (P3.Z - P0.Z) * t;
            return new PrecisedVertex(x, y, z);
        } 
        
        private PrecisedVertex CrossPointVert(double z)
        {
            PrecisedVertex endVector, beginVector;
            PrecisedVertex[] P1, P2;
            GetControlvectors(out P1, out P2, out beginVector, out endVector);
            double a, b;
            double eps = 0.0001;
            int i = 0; double root = -999, X, Y, Z;
            if (z < vertices[0].Z + eps)
            {

                root = (z - vertices[0].Z) / beginVector.Z;
                X = vertices[0].X + beginVector.X * root;
                Y = vertices[0].Y + beginVector.Y * root;
                return new PrecisedVertex(X, Y, z);

            }
            else
                if (z > vertices[4].Z - eps)
                {
                    root = (z - vertices[4].Z) / endVector.Z;
                    X = vertices[4].X + endVector.X * root;
                    Y = vertices[4].Y + endVector.Y * root;
                    return new PrecisedVertex(X, Y, z);

                }
                else if (z < vertices[2].Z)
                {
                    i = (z < vertices[1].Z) ? 0 : 1;
                    root = (z - vertices[i].Z) / (vertices[i + 1].Z - vertices[i].Z);
                }
                else
                {
                    i = (z < vertices[3].Z) ? 2 : 3;
                    root = (z - vertices[i].Z) / (vertices[i + 1].Z - vertices[i].Z);
                }

            if (root == -999) return new PrecisedVertex(-999, -999, -999);
            return BezierCurve(vertices[i], P1[i], P2[i], vertices[i + 1], root);

        }
        /// <summary>
        /// возвращает точку на параметризованной кривой
        ///параметр от 0 до 1. Кривая параметризируется в 
        /// соответствии с координатой z.
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        public PrecisedVertex ParametricByZ(double z)
        {
            return CrossPointVert(z);
        }
        /// <summary>
        /// возвращает точку на параметризованной кривой
        /// параметр от 0 до 1. Кривая параметризируется в 
        /// соответствии с координатой z
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public PrecisedVertex ParametricByT(double t) 
        {
             PrecisedVertex endVector, beginVector;
             PrecisedVertex[] P1, P2;
             GetControlvectors(out P1, out P2, out beginVector, out endVector);

             double[] res = new double[3];
             int i = 0;
             if (t == -999) return new PrecisedVertex(-999, -999, -999);
             if (t <= 0)
             {

                 res[0] = vertices[0].X + beginVector.X * 4 * t;
                 res[1] = vertices[0].Y + beginVector.Y * 4 * t;
                 res[2] = vertices[0].Z + beginVector.Z * 4 * t;
                 return new PrecisedVertex(res[0], res[1], res[2]);
             }
             else if (t >= 1)
             {
                 res[0] = vertices[4].X + endVector.X * 4 * (t - 1);
                 res[1] = vertices[4].Y + endVector.Y * 4 * (t - 1);
                 res[2] = vertices[4].Z + endVector.Z * 4 * (t - 1);
                 return new PrecisedVertex(res[0], res[1], res[2]);
             }
             else if (t < 0.5)
             {
                 i = (t < 0.25) ? 0 : 1;
             }
             else
             {
                 i = (t < 0.75) ? 2 : 3;
             }


             return BezierCurve(vertices[i], P1[i], P2[i], vertices[i + 1], t * 4 - i);
        }
        
        private static double PointOnSurface(BaseSurface surface2, double x, double y)
        {
            double z = -999;

            int i0 = (int)Math.Floor((x - surface2.Xmin) / surface2.XStep);
            int i1 = i0 + 1;
            // int ix = (Math.Abs(x - xmin2 - dx2 * c) < Math.Abs(x - xmin2 - dx2 * f)) ? c : f;

            if (i1 < 0 || i0 > (surface2.Nx - 1))
            {
                return z;
            }



            int j0 = (int)Math.Floor((y - surface2.Ymin) / surface2.YStep);
            int j1 = j0 + 1;
            // int iy= (Math.Abs(y - ymin2 - dy2 * c) < Math.Abs(y - ymin2 - dy2 * f)) ? c : f;

            if (j1 < 0 || j0 > (surface2.Ny - 1)) { return z; }

            if (i0 < 0) i0 = i1;
            if (i1 > (surface2.Nx - 1)) i1 = i0;
            if (j0 < 0) j0 = j1;
            if (j1 > (surface2.Ny - 1)) j1 = j0;

            if (surface2[i0, j0].Z == surface2.NullValue || surface2[i1, j0].Z == surface2.NullValue || surface2[i1, j1].Z == surface2.NullValue || surface2[i0, j1].Z == surface2.NullValue)
                return z;
            double xx = (x - surface2.Xmin - i0 * surface2.XStep) / surface2.XStep;
            double yy = (y - surface2.Ymin - j0 * surface2.YStep) / surface2.YStep;
            double nx, ny, nz;
            if (yy - xx > 0)
            {

                nx = (-1) * (surface2[i1, j1].Z - surface2[i0, j1].Z);
                ny = (surface2[i0, j0].Z - surface2[i0, j1].Z) * (1);
                nz = -(-1) * (1);

                z = (nx * (0 - xx) + ny * (1 - yy)) / nz + surface2[i0, j1].Z;
            }
            else
            {
                nx = -(1) * (surface2[i0, j0].Z - surface2[i1, j0].Z);
                ny = -(surface2[i1, j1].Z - surface2[i1, j0].Z) * (-1);
                nz = (1) * (-1);

                z = (nx * (1 - xx) + ny * (0 - yy)) / nz + surface2[i1, j0].Z;
            }

            return z;

        }
        /// <summary>
        /// возвращает параметр t, при котором пилар пересекает плоскость z
        /// </summary>
        /// <param name="z"></param>
        /// <returns></returns>
        private double CrossPoint(double z)
        {
            PrecisedVertex endVector, beginVector;
            PrecisedVertex[] P1, P2;
            GetControlvectors(out P1, out P2, out beginVector, out endVector);
            double a, b;
            double eps = 0.0001;
            int i = 0;
            if (z < vertices[0].Z + eps)
            {


                return 0.25 * (z - vertices[0].Z) / beginVector.Z;

            }
            else
                if (z > vertices[4].Z - eps)
                {
                    return 1 + 0.25 * (z - vertices[4].Z) / endVector.Z; ;

                }
                else if (z < vertices[2].Z)
                {
                    i = (z < vertices[1].Z) ? 0 : 1;
                }
                else
                {
                    i = (z < vertices[3].Z) ? 2 : 3;
                }

            double root = (z - vertices[i].Z) / (vertices[i + 1].Z - vertices[i].Z);

            if (root != -999) root = root * 0.25 + 0.25 * i;
            return root;
        }
        
        public PrecisedVertex CrossSurfaceExact (BaseSurface surface)
        {

            double t0 = CrossPoint(surface.Min.Z);
            double t1 = CrossPoint(surface.Max.Z);
            int tau0=Math.Max(0,(int)Math.Floor(t0*4));
            int tau1=Math.Min(4,(int)Math.Ceiling(t1*4));
            PrecisedVertex res=new PrecisedVertex(-999,-999,-999);
            for ( int t = tau0; t < tau1;t++ )
                if (CrossSurAndLine(surface, vertices[t], vertices[t+1],out  res))break;
         
            return res;
        }
        
        bool CrossSurAndLine(BaseSurface surface, PrecisedVertex A, PrecisedVertex B, out PrecisedVertex res)
        { 
            
            res = new PrecisedVertex(-999, -999, -999);
            List<double> Intersect = new List<double>();
            List<sbyte> greatThanZ=  new List<sbyte>();
            List<Index3D> surfInd=new List<Index3D>();
            List<int> Ind=new List<int>();
            int nx = surface.Nx;
            int ny = surface.Ny;
            int imin, jmin, imax, jmax,i0,j0;
            double x0, y0,z0,z, a, t, t0,b,c,d;
            imin = (int)Math.Floor((A.X - surface.Xmin) / surface.XStep);
            imax = (int)Math.Floor((B.X - surface.Xmin) / surface.XStep);
            jmin = (int)Math.Floor((A.Y - surface.Ymin) / surface.YStep);
            jmax = (int)Math.Floor((B.Y - surface.Ymin) / surface.YStep);
            if (jmax >= 0 && jmax < ny-1 && imax >= 0 && imax < nx-1)
            {
                Ind.Add(greatThanZ.Count);
                Intersect.Add(1);
                surfInd.Add(new Index3D(imax, jmax, 1));
                greatThanZ.Add(1);
            }
            if (jmin >= 0 && jmin < ny - 1 && imin >= 0 && imin < nx - 1)
            {
                Ind.Add(greatThanZ.Count);
                Intersect.Add(0);
                surfInd.Add(new Index3D(imin, jmin, 1));
                greatThanZ.Add(1);
            }
            a = B.X - A.X;
            if (imin > imax) { i0 = imin; imin = imax; imax = i0; }
            imin = Math.Max(0,imin);
            imax = Math.Min(nx - 1, imax);
            if ( a != 0)
            {

                for (int l = imin + 1; l <= imax; l++)
                {
                    x0 = l * surface.XStep + surface.Xmin;
                    t = (x0 - A.X) / a;
                    if (t >= 0 && t <= 1)
                    {
                        y0 = A.Y + (B.Y - A.Y) * t;
                        if (y0 >= surface.Min.Y && y0 <= surface.Max.Y)
                        {
                            z0 = A.Z + (B.Z - A.Z) * t;
                            i0 = l;
                            y0 = (y0 - surface.Ymin) / surface.YStep;
                            j0 = (int)Math.Floor(y0);
                            Intersect.Add(t);
                            t = (y0 - j0);
                            Ind.Add(greatThanZ.Count);
                          
                             if (j0 == surface.Ny - 1) { j0--; t = 1; }
                            
                             if (surface[i0, j0].Z == surface.NullValue || surface[i0, j0 + 1].Z == surface.NullValue)
                                 greatThanZ.Add(0);
                             else
                             {
                                 z = surface[i0, j0].Z + t * (surface[i0 , j0+1].Z - surface[i0, j0].Z) - z0;
                                 if (z == 0)
                                 {
                                     t = Intersect[Intersect.Count - 1];
                                     res = new PrecisedVertex(A.X + (B.X - A.X) * t, A.Y + (B.Y - A.Y) * t, A.Z + (B.Z - A.Z) * t);
                                     return true;
                                 }
                                 if (z > 0) greatThanZ.Add(1);
                                 else greatThanZ.Add(-1);
                             }
                            
                             if (i0 > 0)
                             { 
                                 if (i0<nx-1)
                                 surfInd.Add(new Index3D(i0, j0,2));
                                 else surfInd.Add(new Index3D(i0, j0, 0));
                               
                             }
                             else surfInd.Add(new Index3D(i0, j0, 1));
                        }
                    }
                }
            }

            a = B.Y - A.Y;
            if (jmin > jmax) { j0 = jmin; jmin = jmax; jmax = j0; }
            jmin = Math.Max(0, jmin);
           jmax = Math.Min(ny - 1, jmax);
            if (a!=0)
            {
                for (int l = jmin + 1; l <= jmax; l++)
                {
                    y0 = l * surface.YStep + surface.Ymin;
                    t = (y0 - A.Y) / a;
                    if (t >= 0 && t <= 1)
                    {
                        x0 = A.X + (B.X - A.X) * t;
                        if (x0 >= surface.Min.X && x0 <= surface.Max.X)
                        {
                            z0 = A.Z + (B.Z - A.Z) * t;
                            j0 = l;
                            x0 = (x0 - surface.Xmin) / surface.XStep;
                            i0 = (int)Math.Floor(x0);
                            Intersect.Add(t);
                            t = (x0 - i0);
                            Ind.Add(greatThanZ.Count);
                             if (i0 == surface.Nx - 1) { i0--; t = 1; }
                             if (surface[i0, j0].Z == surface.NullValue || surface[i0+1, j0].Z == surface.NullValue)
                                 greatThanZ.Add(0);
                             else
                             {
                               z= surface[i0, j0].Z + t * (surface[i0 + 1, j0].Z - surface[i0, j0].Z)-z0;
                                 if (z==0)
                                 {
                                     t=Intersect[Intersect.Count-1];
                                     res=new PrecisedVertex( A.X + (B.X - A.X) *t, A.Y + (B.Y - A.Y) *t, A.Z + (B.Z - A.Z) *t);
                                     return true;
                                 }
                                    if(z>0) greatThanZ.Add(1);
                                 else greatThanZ.Add(-1);
                             }
                             if (j0 > 0)
                             {
                                 if (j0 < ny - 1)
                                     surfInd.Add(new Index3D(i0, j0, 5));
                                 else surfInd.Add(new Index3D(i0, j0, 3));

                             }
                             else surfInd.Add(new Index3D(i0, j0,4));
                          

                        }
                    }
                }
            }
            if (Ind.Count == 0) return false;
            double[] intersect = Intersect.ToArray();
            Intersect.Clear();
            Vertex A0,B0,C0;
            int[] ind = Ind.ToArray();
            Ind.Clear();
            Array.Sort(intersect, ind);
            double v=0 ;
            int iprev1=-1, jprev1=-1;
            int iprev2 = -1, jprev2 = -1;
            for ( int i = 0; i < ind.Length; i++)
            {

                int k=ind[i];
                i0=surfInd[k].I;
                j0=surfInd[k].J;

                if (intersect[i] == 0 || intersect[i] == 1||(i<ind.Length-1&& greatThanZ[k] * greatThanZ[ind[i + 1]] < 0))
                {
                    if (surfInd[k].K % 3 != 0)
                    if (iprev1 != i0 || jprev1 != j0)
                    {
                        A0 = surface[i0, j0];
                        C0 = surface[i0 + 1, j0 + 1];
                        if (A0.Z != surface.NullValue && C0.Z != surface.NullValue)
                        {


                            B0 = surface[i0 + 1, j0];
                            if (B0.Z != surface.NullValue)
                            {
                                TimeZYX.Model.Utils.WellTracerSDB.PlaneCoeff(1, A0, B0, C0, out a, out b, out c, out d);
                                t = TimeZYX.Model.Utils.WellTracerSDB.PlaneIntersect(a, b, c, d, A, B);
                                if (t >= 0 && t <= 1)
                                {
                                    res = new PrecisedVertex(A.X + t * (B.X - A.X), A.Y + t * (B.Y - A.Y), A.Z + t * (B.Z - A.Z));
                                    if (TimeZYX.Model.Utils.WellTracerSDB.inPlane(A0, B0, C0, res)) return true;
                                }

                            }



                            B0 = surface[i0, j0 + 1];
                            if (B0.Z != surface.NullValue)
                            {
                                TimeZYX.Model.Utils.WellTracerSDB.PlaneCoeff(1, A0, B0, C0, out a, out b, out c, out d);
                                t = TimeZYX.Model.Utils.WellTracerSDB.PlaneIntersect(a, b, c, d, A, B);
                                if (t >= 0 && t <= 1)
                                {
                                    res = new PrecisedVertex(A.X + t * (B.X - A.X), A.Y + t * (B.Y - A.Y), A.Z + t * (B.Z - A.Z));
                                    if (TimeZYX.Model.Utils.WellTracerSDB.inPlane(A0, B0, C0, res)) return true;
                                }
                            }

                        }
                        iprev1 = i0; jprev1 = j0;

                    }
                    if (surfInd[k].K % 3 != 1)
                    {
                        if (surfInd[k].K<3)i0--;else j0--;
                        if (iprev2 != i0 || jprev2 != j0)
                        {
                            A0 = surface[i0, j0];
                            C0 = surface[i0 + 1, j0 + 1];
                            if (A0.Z != surface.NullValue && C0.Z != surface.NullValue)
                            {


                                B0 = surface[i0 + 1, j0];
                                if (B0.Z != surface.NullValue)
                                {
                                    TimeZYX.Model.Utils.WellTracerSDB.PlaneCoeff(1, A0, B0, C0, out a, out b, out c, out d);
                                    t = TimeZYX.Model.Utils.WellTracerSDB.PlaneIntersect(a, b, c, d, A, B);
                                    if (t >= 0 && t <= 1)
                                    {
                                        res = new PrecisedVertex(A.X + t * (B.X - A.X), A.Y + t * (B.Y - A.Y), A.Z + t * (B.Z - A.Z));
                                        if (TimeZYX.Model.Utils.WellTracerSDB.inPlane(A0, B0, C0, res)) return true;
                                    }

                                }



                                B0 = surface[i0, j0 + 1];
                                if (B0.Z != surface.NullValue)
                                {
                                    TimeZYX.Model.Utils.WellTracerSDB.PlaneCoeff(1, A0, B0, C0, out a, out b, out c, out d);
                                    t = TimeZYX.Model.Utils.WellTracerSDB.PlaneIntersect(a, b, c, d, A, B);
                                    if (t >= 0 && t <= 1)
                                    {
                                        res = new PrecisedVertex(A.X + t * (B.X - A.X), A.Y + t * (B.Y - A.Y), A.Z + t * (B.Z - A.Z));
                                        if (TimeZYX.Model.Utils.WellTracerSDB.inPlane(A0, B0, C0, res)) return true;
                                    }
                                }

                            }
                            iprev2 = i0; jprev2 = j0;
                        }
                    }


                     
                }
            
            }
            res = new PrecisedVertex(-999, -999, -999);
            return false;
        }
        
        public double CrossSurface(BaseSurface surface)
        {

            double t0 = CrossPoint(surface.Min.Z);
            double t1 = CrossPoint(surface.Max.Z);
            if (t0 == -999 || t1 == -999) return -999;
            return CrossSur(surface, t0, t1);
        }
        
        public void SetImpl(PillarImpl pillarImple)
        {
            impl = pillarImple;
            impl.ConvertPillar(vertices);
        }

        public int Count { get{return impl.Count;} }

        public PrecisedVertex this[int index]
        {
            get { return GetVertex(index); }
            set { SetVertex(value, index); }
        }

        #region IMinMax<PrecisedVertex> Members

        public PrecisedVertex Min
        {
            get { return currentMin; }
        }

        public PrecisedVertex Max
        {
            get { return currentMax; }
        }

        #endregion

        #region IMoveable Members

        public void Move(double dx, double dy, double dz)
        {
            throw new NotImplementedException();
        }
        public bool CanMove(double dx, double dy, double dz)
        {
            return true;
        }
        #endregion

        #region ICloneable Members

        public Pillar Copy()
        {
            return (Pillar)Clone();
        }

        public override object Clone()
        {
            PrecisedVertex[] vertices = new PrecisedVertex[this.vertices.Count];

            this.vertices.CopyTo(vertices,0);
            return new Pillar(Project, Title, vertices, (PillarImpl)(Activator.CreateInstance(impl.GetType())));
        }

        #endregion
    }
}
