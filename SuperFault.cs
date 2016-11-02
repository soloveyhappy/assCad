using System;
using System.Collections.Generic;
using TimeZYX.Core.Base;
using TimeZYX.Model.Spatial;
using TimeZYX.Core.Attributes;
using TimeZYX.Model.Properties;
using System.Drawing;
using System.Collections.ObjectModel;
using TimeZYX.Core.Tools;
using TimeZYX.Model.Utils;
using TimeZYX.Model.Geology.Pillars;

namespace TimeZYX.Model.Geology.FaultsModel
{

    public enum BoundaryFault
    {
        Constant,
        Surface
    }

    public class FaultSettings
    {

        public BoundaryFault UpperBoundaryFault { get; set; }
        public BoundaryFault LowerBoundaryFault { get; set; }
        public BaseSurface SurfaceTop { get; set; }
        public BaseSurface SurfaceDown { get; set; }
        public bool StretchToTheBordersOfTheInterval { get; set; }
        public bool CropIsOutOfRange { get; set; }
        public double Proportion { get; set; }
        public double PillarDistance{ get; set; }
        public double FaultHeight { get; set; }
        public double FaultDown { get; set; }
        public double FaultTop { get; set; }
        public bool UseInterval { get; set; }
        public PillarImpl Impl { get; set; }
        public SuperFault.FaultType FaultType { get; set; }

        public FaultSettings()
        {
            PillarDistance = 300;
            Proportion = 100;
            FaultTop = 1700;
            FaultDown = 2000;
            UseInterval = true;
            FaultHeight = 300;
            Impl = new ListricPillarImpl();
            FaultType = SuperFault.FaultType.sbros;
        }
    }

    [Serializable]
    [DescriptionLocalized("Fault", typeof(Resources))]
    [Image(typeof(Resources), "FaultImage", "FaultImage")]
    public class SuperFault : DataList<Pillar>, IMoveable, IMinMax<PrecisedVertex>, ICompositor<SuperFault>
    {
        private static readonly Image nodeImage;

        public override Image Image
        {
            get { return nodeImage; }
        }

        static SuperFault()
        {
            nodeImage = Resources.FaultImage;
        }

        public SuperFault(Project prj, string title)
            : base(prj, title)
        {
            Init();
        }

        public SuperFault(Project prj, string title, int cap)
            : base(prj, title, cap)
        {
            Init();
        }
        
        public SuperFault(Project prj, string title, IEnumerable<Pillar> pillars)
            : base(prj, title, pillars)
        {
            Init();
            RecalcMinMax();
        }

        public SuperFault(Project prj, string title, IEnumerable<Pillar> pillars, 
                          IEnumerable<Pillar> leftConnPillars, IEnumerable<Pillar> rightConnPillars)
            : this(prj, title, pillars)
        {
            this.leftConnectedPillar = new List<Pillar>(leftConnPillars);
            leftConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(leftConnectedPillar);
            this.rightConnectedPillar = new List<Pillar>(rightConnPillars);
            rightConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(rightConnectedPillar);
        }

        private void RecalcMinMax()
        {
            double xmin = double.MaxValue, ymin = double.MaxValue, zmin = double.MaxValue;
            double xmax = double.MinValue, ymax = double.MinValue, zmax = double.MinValue;

            foreach (Pillar item in items)
            {
                if (xmin > item.Min.X)
                    xmin = item.Min.X;
                if (ymin > item.Min.Y)
                    ymin = item.Min.Y;
                if (zmin > item.Min.Z)
                    zmin = item.Min.Z;

                if (xmax < item.Max.X)
                    xmax = item.Max.X;
                if (ymax < item.Max.Y)
                    ymax = item.Max.Y;
                if (zmax < item.Max.Z)
                    zmax = item.Max.Z;

            }

            currentMin = new PrecisedVertex(xmin, ymin, zmin);
            currentMax = new PrecisedVertex(xmax, ymax, zmax);

        }

       
        private void Init()
        {
            leftConnectedPillar = new List<Pillar>();
            leftConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(leftConnectedPillar);
            
            rightConnectedPillar = new List<Pillar>();
            rightConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(rightConnectedPillar);
        }

        #region IMoveable Members

        public void Move(double dx, double dy, double dz)
        {
            foreach (var item in items)
            {
                item.Move(dx, dy, dz);
            }
        }
        public bool CanMove(double dx, double dy, double dz)
        {
            return true;
        }
        #endregion

        #region IMinMax<PrecisedVertex> Members

        public PrecisedVertex Min
        {
            get 
            {
                return currentMin;
            }
        }

        public PrecisedVertex Max
        {
            get 
            {
                return currentMax;
            }
        }

        #endregion

        /// <summary>
        /// Получить все пиллары с учетом присоединенных
        /// </summary>
        /// <returns>Набор пилларов</returns>
        public IList<Pillar> GetPillarsWithConnected()
        {
            List<Pillar> pillars = new List<Pillar>(Count);
            foreach (var pillar in leftConnectedPillar)
                pillars.Add(pillar);
            foreach (var pillar in items)
                pillars.Add(pillar);
            foreach (var pillar in rightConnectedPillar)
                pillars.Add(pillar);
            return pillars;
        }

        public Pillar LeftPillar
        {
            get
            {
                if (items.Count > 0)
                {
                    return items[0];
                }
                else
                    return null;
            }
        }
        
        public Pillar RightPillar
        {
            get
            {
                if (items.Count > 0)
                {
                    return items[items.Count - 1];
                }
                else
                    return null;
            }
        }

        public ReadOnlyCollection<Pillar> LeftConnectedPillar
        {
            get
            {
                return leftConnectedPillarWrapper;
            }
        }

        public ReadOnlyCollection<Pillar> RightConnectedPillar
        {
            get
            {
                return rightConnectedPillarWrapper;
            }
        }

        /// <summary>
        /// Возвращает пару пиларов first = первый пиллар second = последний пиллар.
        /// </summary>
        public Pair<Pillar, Pillar> EndPillars
        {
            get 
            { 
                Pair<Pillar, Pillar> p = new Pair<Pillar,Pillar>();

                if(items.Count>0)
                {
                    p = new Pair<Pillar, Pillar>(items[0], items[items.Count-1]);
                }

                return p; 
            }
        }

        /// <summary>
        /// Функция возвращает true если переданного пилара нету в коллекции уже соединенных пиларов.
        /// </summary>
        /// <param name="pillar">Проверяемый пиллар.</param>
        /// <returns></returns>
        public bool CanBeConnectedEndPillar(Pillar childPillar)
        {

            // если пилар не содержится в этом разломе
            if (!IsMyPillar(childPillar))
                return false;
            
            // если этот пилар уже с кем-то соединен
            foreach (var item in leftConnectedPillar)
            {
                if (item == childPillar)
                    return false;
            }
            
            foreach (var item in rightConnectedPillar)
            {
                if (item == childPillar)
                    return false;
            }

            return true;
        }

        /// <summary>
        /// Функция возвращает true если переданного пилара нету в коллекции уже соединенных пиларов.
        /// </summary>
        /// <param name="pillar">Проверяемый пиллар.</param>
        /// <returns></returns>
        public bool CanBeConnectedWithCreateNewFault(Pillar childPillar)
        {
            // если пилар содержится в этом разломе
            if (IsMyPillar(childPillar))
                return true;

            return false;
        }



        /// <summary>
        /// Соединяет крайний пилар с sourcePillar пиларом разлома fault.
        /// </summary>
        /// <param name="fault">Разлом содержащий sourcePillar.</param>
        /// <param name="sourcePillar">Пиллар с которым производится соединение.</param>
        /// <param name="childPillar">Крайний пиллар в текущем разломе.</param>
        /// <returns></returns>
        public bool ConnectEndPillar(SuperFault fault, Pillar sourcePillar, Pillar childPillar)
        {
            bool res = false;

            // проверка на то что пилар от этого разлома
            if(!IsMyPillar(childPillar))
                return false;

            if(childPillar.Equals(LeftPillar))
            {
                sourcePillar.Attached = true;
                leftConnectedPillar.Add(sourcePillar);
                return true;
            }
            else if (childPillar.Equals(RightPillar))
            {
                sourcePillar.Attached = true;
                rightConnectedPillar.Add(sourcePillar);
                return true;
            }

            return res;
        }

        /// <summary>
        /// Соединяет пилар с sourcePillar пиларом разлома fault 
        /// с созданием нового пустого разлома с двумя присоединенными пиларами.
        /// </summary>
        /// <param name="fault">Разлом содержащий sourcePillar.</param>
        /// <param name="sourcePillar">Пиллар с которым производится соединение.</param>
        /// <param name="childPillar">Крайний пиллар в текущем разломе.</param>
        /// <returns></returns>
        public SuperFault ConnectWithCreateNewFault(SuperFault fault, Pillar sourcePillar, Pillar childPillar)
        {
            SuperFault res = new SuperFault(this.Project, "fault");
            sourcePillar.Attached = true;
            childPillar.Attached = true;
            res.leftConnectedPillar.Add(sourcePillar);
            res.rightConnectedPillar.Add(childPillar);

            return res;
        }

        private bool IsMyPillar(Pillar p)
        {
            foreach (var item in items)
            {
                if (item.Equals(p))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Возвращает совакупность 
        /// </summary>
        /// <returns></returns>
        public List<Pillar> GetAllPillars()
        { 
            List<Pillar> l = new List<Pillar>();

            l.AddRange(leftConnectedPillar);
            l.AddRange(items);
            l.AddRange(rightConnectedPillar);
            return l;
        }

        #region Перегруженные Функции по изменению списка пиларов (рачет мин и макс)
        public override void Add(Pillar item)
        {
            base.Add(item);
            RecalcMinMax();
        }

        public override void AddRange(IEnumerable<Pillar> collection)
        {
            base.AddRange(collection);
            RecalcMinMax();
        }

        public override void Insert(int index, Pillar item)
        {
            base.Insert(index, item);
            RecalcMinMax();
        }

        public override bool Remove(Pillar item)
        {
            bool res = base.Remove(item);
            RecalcMinMax();
            return res;
        }
        public override void RemoveAt(int index)
        {
            base.RemoveAt(index);
            RecalcMinMax();
        }

        #endregion

        PrecisedVertex[] GetTop()
        {
            PrecisedVertex[] top = new PrecisedVertex[Count];
            int iter = 0;
            foreach (var item in Items)
                top[iter++] = item[0];
            return top;
        }

        PrecisedVertex[] GetBot()
        {
            PrecisedVertex[] bot = new PrecisedVertex[Count];
            int iter = 0;
            foreach (var item in Items)
                bot[iter++] = item[item.Count-1];
            return bot;
        }

        PrecisedVertex[] GetMid()
        {
            PrecisedVertex[] mid = new PrecisedVertex[Count];
            int iter = 0;
            foreach (var item in Items)
                mid[iter++] = item.Center();
            return mid;
        }

        public Pillar AddPillarToEnd(List<int> ids,FaultSettings fsettings, double distance)
        {
            if (ids.Count == 0)
                return null;
            bool isFirst=false;
            ids.Sort();
            if (ids[0] == 0)
                isFirst = true;
            else
            {
                if (ids[ids.Count - 1] != Count - 1)
                    return null;
            }
            List<PrecisedVertex> v = new List<PrecisedVertex>();
            foreach (var item in Items)
            {
                v.Add(item[0]);
            }
            //foreach (var id in ids)
            //{
            //    v.Add(Items[id][0]);
            //}
            PrecisedVertex[] data;
            PrecisedVertex[] P1=new PrecisedVertex[v.Count];
            PrecisedVertex[] P2=new PrecisedVertex[v.Count];
            data = v.ToArray();
            try
            {
                MathFunctions.GetControlPoints(v.ToArray(), out P1, out P2);
                PrecisedVertex top = AddNewPillarToEnd(distance, isFirst, ref data, ref P1, ref P2);
                v.Clear();
                //foreach (var id in ids)
                //{
                //    v.Add(Items[id][Items[id].Count - 1]);
                //}
                foreach (var item in Items)
                {
                    v.Add(item[item.Count - 1]);
                }
                data = v.ToArray();
                MathFunctions.GetControlPoints(v.ToArray(), out P1, out P2);
                PrecisedVertex bot = AddNewPillarToEnd(distance, isFirst, ref data, ref P1, ref P2);
                double len = top.DistanceTo(bot);
                double step = len / 5;
                double vx = (bot.X - top.X) / 4;
                double vy = (bot.Y - top.Y) / 4;
                double vz = (bot.Z - top.Z) / 4;
                PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
                 PrecisedVertex nv = top;
                 List<PrecisedVertex> nlist = new List<PrecisedVertex>(5);
                 nlist.Add(top);
                for (int i = 0; i < 3; i++)
                {
                    nv = nv + vec;
                    nlist.Add(nv);
                }
                nlist.Add(bot);
                Pillar pillar = new Pillar(Project, "Adding", nlist, fsettings.Impl);
                if (isFirst)
                    Insert(0, pillar);
                else
                    Add(pillar);
                //Insert(isFirst ? 0 : Count - 1, pillar);
                return pillar;
            }
            catch
            {
                return null;
            }
            
        }

        public Pillar AddPillarBetween(List<int> ids, FaultSettings fsettings)
        {
            if (ids.Count < 2)
                return null;
            ids.Sort();
            if (ids[1] != ids[0] + 1)
                return null;
            PrecisedVertex top = (Items[ids[0]][0] + Items[ids[1]][0]) / 2.0;
            PrecisedVertex bot = (Items[ids[0]][Items[ids[0]].Count - 1] + Items[ids[1]][Items[ids[1]].Count - 1]) / 2.0;
            double len = top.DistanceTo(bot);
            double step = len / 5;
            double vx = (bot.X - top.X) / 4;
            double vy = (bot.Y - top.Y) / 4;
            double vz = (bot.Z - top.Z) / 4;
            PrecisedVertex vec = new PrecisedVertex(vx, vy, vz);
            PrecisedVertex nv = top;
            List<PrecisedVertex> nlist = new List<PrecisedVertex>(5);
            nlist.Add(top);
            for (int i = 0; i < 3; i++)
            {
                nv = nv + vec;
                nlist.Add(nv);
            }
            nlist.Add(bot);
            Pillar pillar = new Pillar(Project, "Adding", nlist, fsettings.Impl);
            if (pillar != null)
                Insert(ids[1], pillar);
            return pillar;
            //return null;
        }

        static void AddNewPillarBetween(double Inc, int index, ref PrecisedVertex[] data, ref PrecisedVertex[] P1, ref PrecisedVertex[] P2)
        {
            double derX = 0.125 * (-data[index].X - P1[index].X + P2[index].X + data[index + 1].X);
            double derY = 0.125 * (-data[index].Y - P1[index].Y + P2[index].Y + data[index + 1].Y);
            double derZ = 0.125 * (-data[index].Z - P1[index].Z + P2[index].Z + data[index + 1].Z);

            PrecisedVertex V = MathFunctions.BezierCurve(data[index], P1[index], P2[index], data[index + 1], 0.5);

            List<PrecisedVertex> Buf = new List<PrecisedVertex>(data);
            Buf.Insert(index + 1, V);
            data = Buf.ToArray();
            Buf = new List<PrecisedVertex>(P1);
            Buf.Insert(index + 1, new PrecisedVertex());
            P1 = Buf.ToArray();
            Buf = new List<PrecisedVertex>(P2);
            Buf.Insert(index + 1, new PrecisedVertex());
            P2 = Buf.ToArray();
            P2[index + 1] = new PrecisedVertex(P2[index].X, P2[index].Y, P2[index].Z);
            //P1[index + 1].X = (float)(data[index + 1].X + derX);
            //P1[index + 1].Y = (float)(data[index + 1].Y + derY);
            //P1[index + 1].Z = (float)(data[index + 1].Z + derZ);
            P1[index + 1] = P1[index + 1].Move((data[index + 1].X + derX),
                                               (data[index + 1].Y + derY), 
                                               (data[index + 1].Z + derZ));
            //P2[index].X = (float)(data[index + 1].X - derX);
            //P2[index].Y = (float)(data[index + 1].Y - derY);
            //P2[index].Z = (float)(data[index + 1].Z - derZ);
            P2[index] = P2[index].Move((data[index + 1].X - derX), (data[index + 1].Y - derY), (data[index + 1].Z - derZ));
            //P2[index + 1].X = data[index + 2].X + (P2[index + 1].X - data[index + 2].X) / 2;
            //P2[index + 1].Y = data[index + 2].Y + (P2[index + 1].Y - data[index + 2].Y) / 2;
            //P2[index + 1].Z = data[index + 2].Z + (P2[index + 1].Z - data[index + 2].Z) / 2;
            P2[index + 1] = P2[index + 1].Move(data[index + 2].X + (P2[index + 1].X - data[index + 2].X) / 2,
                                               data[index + 2].Y + (P2[index + 1].Y - data[index + 2].Y) / 2,
                                               data[index + 2].Z + (P2[index + 1].Z - data[index + 2].Z) / 2);

            //P1[index].X = data[index].X + (P1[index].X - data[index].X) / 2;
            //P1[index].Y = data[index].Y + (P1[index].Y - data[index].Y) / 2;
            //P1[index].Z = data[index].Z + (P1[index].Z - data[index].Z) / 2;
            P1[index] = P1[index].Move(data[index].X + (P1[index].X - data[index].X) / 2,
                                       data[index].Y + (P1[index].Y - data[index].Y) / 2,
                                       data[index].Z + (P1[index].Z - data[index].Z) / 2);
        }
        

        static PrecisedVertex AddNewPillarToEnd(double Inc, bool isFirst, ref PrecisedVertex[] data, ref PrecisedVertex[] P1, ref PrecisedVertex[] P2)
        {
            double x, y, z, ax, ay, az, l;

            if (isFirst)
            {
                List<PrecisedVertex> Buf;
                l = Math.Sqrt((data[0].X - P1[0].X) * (data[0].X - P1[0].X) + (data[0].Y - P1[0].Y) * (data[0].Y - P1[0].Y) + (data[0].Z - P1[0].Z) * (data[0].Z - P1[0].Z));
                ax = (data[0].X - P1[0].X) * Inc / l;
                ay = (data[0].Y - P1[0].Y) * Inc / l;
                az = (data[0].Z - P1[0].Z) * Inc / l;
                x = data[0].X + ax / 3;
                y = data[0].Y + ay / 3;
                z = data[0].Z + az / 3;
                Buf = new List<PrecisedVertex>(P2);
                Buf.Insert(0, new PrecisedVertex((float)x, (float)y, (float)z));
                P2 = Buf.ToArray();
                x = data[0].X + 2 * ax / 3;
                y = data[0].Y + 2 * ay / 3;
                z = data[0].Z + 2 * az / 3;
                Buf = new List<PrecisedVertex>(P1);
                Buf.Insert(0, new PrecisedVertex((float)x, (float)y, (float)z));
                P1 = Buf.ToArray();
                x = data[0].X + ax;
                y = data[0].Y + ay;
                z = data[0].Z + az;
                Buf = new List<PrecisedVertex>(data);
                Buf.Insert(0, new PrecisedVertex(x, y, z));
                data = Buf.ToArray();
                return new PrecisedVertex(x, y, z);
            }
            else
            {
                int N = data.Length - 2;
                l = Math.Sqrt((data[N + 1].X - P2[N].X) * (data[N + 1].X - P2[N].X) + (data[N + 1].Y - P2[N].Y) * (data[N + 1].Y - P2[N].Y) + (data[N + 1].Z - P2[N].Z) * (data[N + 1].Z - P2[N].Z));
                ax = (data[N + 1].X - P2[N].X) * Inc / l;
                ay = (data[N + 1].Y - P2[N].Y) * Inc / l;
                az = (data[N + 1].Z - P2[N].Z) * Inc / l;
                x = data[N + 1].X + ax / 3;
                y = data[N + 1].Y + ay / 3;
                z = data[N + 1].Z + az / 3;
                Array.Resize(ref P1, P1.Length + 1);
                P1[N + 1] = new PrecisedVertex((float)x, (float)y, (float)z);
                x = data[N + 1].X + 2 * ax / 3;
                y = data[N + 1].Y + 2 * ay / 3;
                z = data[N + 1].Z + 2 * az / 3;
                Array.Resize(ref P2, P2.Length + 1);
                P2[N + 1] = new PrecisedVertex((float)x, (float)y, (float)z);
                x = data[N + 1].X + ax;
                y = data[N + 1].Y + ay;
                z = data[N + 1].Z + az;

                Array.Resize(ref data, data.Length + 1);
                data[N + 2] = new PrecisedVertex(x, y, z);
                return new PrecisedVertex(x, y, z);
            }


        }

        
        PrecisedVertex currentMin;
        PrecisedVertex currentMax;
        List<Pillar> rightConnectedPillar;
        List<Pillar> leftConnectedPillar;
        ReadOnlyCollection<Pillar> rightConnectedPillarWrapper;
        ReadOnlyCollection<Pillar> leftConnectedPillarWrapper;
        FaultType faultType = FaultType.sbros;

        public enum FaultType
        {
            sbros = 1,
            vbros = -1
        }

        public FaultType FType { get{return faultType;} set{faultType=value;} }

        #region ICompositor<SuperFault> Members

        public virtual DataCompositor<SuperFault> Compositor
        {
            get { return null; }
        }

        #endregion

        public override object Clone()
        {
            SuperFault fault = new SuperFault(this.Project, this.Title);
            fault.items = ((DataList<Pillar>)base.Clone()).Items;
            fault.currentMax = this.currentMax;
            fault.currentMin = this.currentMin;
            fault.faultType = this.faultType;

            List<Pillar> pillars = new List<Pillar>(rightConnectedPillar.Count);
            foreach (Pillar pillar in rightConnectedPillar)
                pillars.Add((Pillar)pillar.Clone());
            fault.rightConnectedPillar = pillars;

            pillars = new List<Pillar>(leftConnectedPillar.Count);
            foreach (Pillar pillar in leftConnectedPillar)
                pillars.Add((Pillar)pillar.Clone());
            fault.leftConnectedPillar = pillars;

            pillars = new List<Pillar>(leftConnectedPillarWrapper.Count);
            foreach (Pillar pillar in leftConnectedPillarWrapper)
                pillars.Add((Pillar)pillar.Clone());
            fault.leftConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(pillars);

            pillars = new List<Pillar>(rightConnectedPillarWrapper.Count);
            foreach (Pillar pillar in rightConnectedPillarWrapper)
                pillars.Add((Pillar)pillar.Clone());
            fault.rightConnectedPillarWrapper = new ReadOnlyCollection<Pillar>(pillars);

            return fault;
        }
    }
}


