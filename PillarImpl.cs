using System;
using System.Collections.Generic;
using TimeZYX.Model.Spatial;

namespace TimeZYX.Model.Geology.Pillars
{
    [Serializable]
    public abstract class PillarImpl
    {
        public abstract PrecisedVertex GetCenter(Pillar pillar);
        public abstract PrecisedVertex GetCenter(IList<PrecisedVertex> vertices);
        public abstract PrecisedVertex GetVertex(IList<PrecisedVertex> vertices, int index);
        public abstract IList<PrecisedVertex> GetVertices(IList<PrecisedVertex> vertices);
        public abstract void ConvertPillar(IList<PrecisedVertex> vertices);
        public virtual int Count { get{return 2;} }
        public abstract PillarKind Kind { get; }

        public abstract void SetVertex(IList<PrecisedVertex> vertices, PrecisedVertex newVertex, int index);
        public abstract void SetVertexes(IList<PrecisedVertex> vertices, PrecisedVertex[] newVertexes);
        /// <summary>
        /// возвращает точку на параметризованной кривой
        /// параметр от 0 до 1. Кривая параметризируется в 
        /// соответствии с координатой z
        /// 
        /// Нужно уточнить нужна ли эта функция вне этого класса?
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public virtual PrecisedVertex ParametricByT(double t)
        {
            return new PrecisedVertex();
        }
        public virtual double CrossPoint(double z)
        {
             return 0;
        }
    }

    [Serializable]
    public class LinearPillarImpl : PillarImpl
    {
        public override PrecisedVertex GetCenter(IList<PrecisedVertex> vertices)
        {
            return (vertices[0] + vertices[vertices.Count - 1]) / 2;
        }

        public override IList<PrecisedVertex> GetVertices(IList<PrecisedVertex> vertices)
        {
            return new List<PrecisedVertex> { vertices[0], vertices[vertices.Count-1]};
        }

        public override void SetVertex(IList<PrecisedVertex> vertices, PrecisedVertex newVertex, int index)
        {
            if (index == 0)
            {
                if (newVertex.Z < vertices[vertices.Count - 1].Z)
                    vertices[0] = newVertex;
            }
            else
            {
                if (newVertex.Z > vertices[0].Z)
                    vertices[vertices.Count - 1] = newVertex;
            }
            double len = vertices[0].DistanceTo(vertices[vertices.Count - 1]);
            double step = len / vertices.Count;
            double vx = (vertices[vertices.Count - 1].X - vertices[0].X) / (vertices.Count - 1);
            double vy = (vertices[vertices.Count - 1].Y - vertices[0].Y) / (vertices.Count - 1);
            double vz = (vertices[vertices.Count - 1].Z - vertices[0].Z) / (vertices.Count - 1);
            PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
            for (int i = 1; i < vertices.Count - 1; i++)
                vertices[i] = vertices[i - 1] + vec;
        }
        public override void SetVertexes(IList<PrecisedVertex> vertices, PrecisedVertex[] newVertexes)
        {
            if (newVertexes.Length != 2) return;
            if (newVertexes[0].Z > newVertexes[1].Z)
            { vertices[0] = newVertexes[1];  vertices[4] = newVertexes[0];}
            else { vertices[0] = newVertexes[0]; vertices[4] = newVertexes[1]; }
           
        }
        public override PrecisedVertex GetCenter(Pillar pillar)
        {
            throw new NotImplementedException();
        }

        public override PrecisedVertex GetVertex(IList<PrecisedVertex> vertices, int index)
        {
            return index == 0 ? vertices[0] : vertices[vertices.Count - 1]; 
        }

        public override void ConvertPillar(IList<PrecisedVertex> vertices)
        {
            double len = vertices[0].DistanceTo(vertices[vertices.Count - 1]);
            double step = len / vertices.Count;
            double vx = (vertices[vertices.Count - 1].X - vertices[0].X) / (vertices.Count-1);
            double vy = (vertices[vertices.Count - 1].Y - vertices[0].Y) / (vertices.Count-1);
            double vz = (vertices[vertices.Count - 1].Z - vertices[0].Z) / (vertices.Count-1);
            PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
            for (int i = 1; i < vertices.Count-1; i++)
            {
                vertices[i] = vertices[i - 1] + vec;
            }
            //buf = new PrecisedVertex[2];
            //vertices[0] = items[0];
            //vertices[vertices.Count-1] = items[items.Count - 1];
        }

        public override int Count { get { return 2; } }

        public override PillarKind Kind { get { return PillarKind.LINEAR; } }
    }

    [Serializable]
    public class ListricPillarImpl : PillarImpl
    {
        public override PrecisedVertex GetCenter(Pillar pillar)
        {
            throw new NotImplementedException();
        }

        public override PrecisedVertex GetCenter(IList<PrecisedVertex> vertices)
        {
            return vertices[2];
        }

        public override PrecisedVertex GetVertex(IList<PrecisedVertex> vertices, int index)
        {
            return index == 0 ? vertices[0] : (index == 1 ? vertices[2] : vertices[4]);

            //if (index == 0)
            //    return vertices[0];
            //else if (index == 1)
            //    return vertices[2];
            //else
            //    return vertices[4];
        }

        public override IList<PrecisedVertex> GetVertices(IList<PrecisedVertex> vertices)
        {
            return new List<PrecisedVertex> { vertices[0], vertices[2], vertices[4] };
        }

        public override void ConvertPillar(IList<PrecisedVertex> vertices)
        {
            //LinearPillarImpl tmp = new LinearPillarImpl();
            //tmp.ConvertPillar(vertices);
           // vertices[1] =new PrecisedVertex();
        }

        public override void SetVertex(IList<PrecisedVertex> vertices, PrecisedVertex newVertex, int index)
        {
            if (index == 0)
            {
                if(newVertex.Z<vertices[2].Z)
                vertices[0] = newVertex;
            }
            else if (index == 1)
            {
                if (newVertex.Z > vertices[0].Z && newVertex.Z < vertices[4].Z)
                vertices[2] = newVertex;
            }
            else
            {
                if (newVertex.Z > vertices[2].Z)
                vertices[4] = newVertex;
            }

            PrecisedVertex[] p1, p2;
            List<PrecisedVertex> buf = new List<PrecisedVertex>(this.GetVertices(vertices));
            if ((buf.Count - 1) <= 0) return;
            TimeZYX.Model.Utils.MathFunctions.GetControlPoints(buf.ToArray(), out p1, out p2);
            //tmp[0] = buf[0];
            //tmp[4] = buf[buf.Count - 1];
            //for (int k = 1; k < 4; k++)
            //{
                double t = 0.25 * 1;
                int p = (int)(t * (buf.Count - 1));
                t = (t - p * 1.0 / (buf.Count - 1));
                vertices[1] = TimeZYX.Model.Utils.MathFunctions.BezierCurve(buf[p], p1[p], p2[p], buf[p + 1], t);

                t = 0.25 * 3;
                p = (int)(t * (buf.Count - 1));
                t = (t - p * 1.0 / (buf.Count - 1));
                vertices[3] = TimeZYX.Model.Utils.MathFunctions.BezierCurve(buf[p], p1[p], p2[p], buf[p + 1], t);

            //}
        }
        public override void SetVertexes(IList<PrecisedVertex> vertices, PrecisedVertex[] newVertexes)
        {
          if (newVertexes.Length!=3) return;
            if (newVertexes[0].Z > vertices[2].Z)
                vertices[0] = new PrecisedVertex(newVertexes[0].X, newVertexes[0].Y, vertices[2].Z);
            else vertices[0] = newVertexes[0];
            if (newVertexes[2].Z < vertices[2].Z)
                vertices[4] = new PrecisedVertex(newVertexes[2].X, newVertexes[2].Y, vertices[2].Z);
            else vertices[4] = newVertexes[2];
            PrecisedVertex[] p1, p2;
            List<PrecisedVertex> buf = new List<PrecisedVertex>(this.GetVertices(vertices));
            if ((buf.Count - 1) <= 0) return;
            TimeZYX.Model.Utils.MathFunctions.GetControlPoints(buf.ToArray(), out p1, out p2);
            //tmp[0] = buf[0];
            //tmp[4] = buf[buf.Count - 1];
            //for (int k = 1; k < 4; k++)
            //{
            double t = 0.25 * 1;
            int p = (int)(t * (buf.Count - 1));
            t = (t - p * 1.0 / (buf.Count - 1));
            vertices[1] = TimeZYX.Model.Utils.MathFunctions.BezierCurve(buf[p], p1[p], p2[p], buf[p + 1], t);

            t = 0.25 * 3;
            p = (int)(t * (buf.Count - 1));
            t = (t - p * 1.0 / (buf.Count - 1));
            vertices[3] = TimeZYX.Model.Utils.MathFunctions.BezierCurve(buf[p], p1[p], p2[p], buf[p + 1], t);

            //}
        }
        public override int Count { get { return 3; } }

        public override PillarKind Kind { get { return PillarKind.LISTRIC; } }
    }

    [Serializable]
    public class VerticalPillarImpl : PillarImpl
    {
        public override PrecisedVertex GetCenter(Pillar pillar)
        {
            throw new NotImplementedException();
        }

        public override PrecisedVertex GetCenter(IList<PrecisedVertex> vertices)
        {
            return (vertices[0] + vertices[vertices.Count - 1]) / 2;
        }

        public override PrecisedVertex GetVertex(IList<PrecisedVertex> vertices, int index)
        {
            return index == 0 ? vertices[0] : vertices[vertices.Count - 1]; 
        }

        public override IList<PrecisedVertex> GetVertices(IList<PrecisedVertex> vertices)
        {
            return new List<PrecisedVertex> { vertices[0], vertices[vertices.Count - 1] };
        }

        public override void ConvertPillar(IList<PrecisedVertex> vertices)
        {
            PrecisedVertex center = GetCenter(vertices);
            vertices[0] = new PrecisedVertex(center.X, center.Y, vertices[0].Z);
            vertices[vertices.Count - 1] = new PrecisedVertex(center.X, center.Y, vertices[vertices.Count - 1].Z);
            double len = vertices[0].DistanceTo(vertices[vertices.Count - 1]);
            double step = len / vertices.Count;
            double vx = (vertices[vertices.Count - 1].X - vertices[0].X) / (vertices.Count-1);
            double vy = (vertices[vertices.Count - 1].Y - vertices[0].Y) / (vertices.Count-1);
            double vz = (vertices[vertices.Count - 1].Z - vertices[0].Z) / (vertices.Count-1);
            PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
            for (int i = 1; i < vertices.Count - 1; i++)
            {
                vertices[i] = vertices[i - 1] + vec;
            }
        }

        public override void SetVertex(IList<PrecisedVertex> vertices, PrecisedVertex newVertex, int index)
        {

            if (index == 0)
            {
                if (newVertex.Z < vertices[vertices.Count - 1].Z)
                    vertices[0] = newVertex;
            }
            else
            {
                if (newVertex.Z > vertices[0].Z)
                vertices[vertices.Count - 1] = newVertex;
            }

            double len = vertices[0].DistanceTo(vertices[vertices.Count - 1]);
            double step = len / vertices.Count;
            double vx = (vertices[vertices.Count - 1].X - vertices[0].X) / (vertices.Count - 1);
            double vy = (vertices[vertices.Count - 1].Y - vertices[0].Y) / (vertices.Count - 1);
            double vz = (vertices[vertices.Count - 1].Z - vertices[0].Z) / (vertices.Count - 1);
            PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
            for (int i = 1; i < vertices.Count - 1; i++)
            {
                vertices[i] = vertices[i - 1] + vec;
            }
            //ConvertPillar(vertices);
            //if (index == 0)
            //    vertices[0] = newVertex;
            //else if (index == 1)
            //    vertices[2] = newVertex;
            //else
            //    vertices[4] = newVertex;
        }
        public override void SetVertexes(IList<PrecisedVertex> vertices, PrecisedVertex[] newVertex)
        { 
        }
        public override int Count { get { return 2; } }

        public override PillarKind Kind { get { return PillarKind.VERTICAL; } }
    }

    [Serializable]
    public class CurvedPillarImpl : PillarImpl
    {
        public override PrecisedVertex GetCenter(Pillar pillar)
        {
            throw new NotImplementedException();
        }

        public override PrecisedVertex GetCenter(IList<PrecisedVertex> vertices)
        {
            return vertices[2];
        }

        public override PrecisedVertex GetVertex(IList<PrecisedVertex> vertices, int index)
        {
            return vertices[index];
        }

        public override IList<PrecisedVertex> GetVertices(IList<PrecisedVertex> vertices)
        {
            return vertices;
        }

        public override void ConvertPillar(IList<PrecisedVertex> vertices)
        {
            //LinearPillarImpl tmp = new LinearPillarImpl();
            //tmp.ConvertPillar(vertices);
        }

        public override void SetVertex(IList<PrecisedVertex> vertices, PrecisedVertex newVertex, int index)
        {
            if(index ==0)
            {
                if(newVertex.Z<vertices[1].Z)
                    vertices[index] = newVertex;
            }
            else if(index==1)
            {
                if (newVertex.Z > vertices[0].Z && newVertex.Z < vertices[2].Z)
                    vertices[index] = newVertex;
            }
            else if(index==2)
            {
                if (newVertex.Z > vertices[1].Z && newVertex.Z < vertices[3].Z)
                    vertices[index] = newVertex;
            }
            else if(index==3)
            {
                if (newVertex.Z > vertices[2].Z && newVertex.Z < vertices[4].Z)
                    vertices[index] = newVertex;
            }
            else
            {
                 if(newVertex.Z>vertices[3].Z)
                    vertices[index] = newVertex;
            }
            
        }
          public override void SetVertexes(IList<PrecisedVertex> vertices, PrecisedVertex[] newVertexes)
        {
            if (newVertexes.Length != 5) return;
            if (newVertexes[1].Z > vertices[2].Z)
                vertices[1] = new PrecisedVertex(newVertexes[1].X, newVertexes[1].Y, vertices[2].Z);
            else vertices[1] = newVertexes[1];
            if (newVertexes[3].Z < vertices[2].Z)
                vertices[3] = new PrecisedVertex(newVertexes[3].X, newVertexes[3].Y, vertices[2].Z);
            else vertices[3] = newVertexes[3];
            if (newVertexes[0].Z > vertices[1].Z)
                vertices[0] = new PrecisedVertex(newVertexes[0].X, newVertexes[0].Y, vertices[1].Z);
            else vertices[0] = newVertexes[0];
            if (newVertexes[4].Z < vertices[3].Z)
                vertices[4] = new PrecisedVertex(newVertexes[4].X, newVertexes[4].Y, vertices[3].Z);
            else vertices[4] = newVertexes[4];
            
        }
         
        public override int Count { get { return 5; } }

        public override PillarKind Kind { get { return PillarKind.CURVED; } }
    }
}
