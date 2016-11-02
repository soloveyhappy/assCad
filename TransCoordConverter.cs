using System;
using System.Reflection;
using System.Collections.Generic;
using TimeZYX.Model.Wells;
using TimeZYX.Model.Properties;
using System.Resources;

namespace TimeZYX.Model.Geology
{
    /// <summary>
    /// Класс для конвертации из различных геодезических систем
    /// </summary>
    public class TransCoordConverter
    {
        #region Public methods

        /// <summary>
        /// Перевод СК-42 в ПЗ-90
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_266", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromCK42ToPZ90(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                             25, -141, -80,
                             0, -0.35, -0.66, 0,
                             6378245, 298.3,         // Эллипсоид Красовского
                             6378136, 298.257839303, // Эллипсоид ПЗ-90
                             direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378136, 298.257839303, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        /// <summary>
        /// Перевод СК-42 в ПЗ-90.02
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_267", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromCK42ToPZ9002(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                             23.93, -141.03, -79.98,
                             0, 0, 0, -0.22E-06,
                             6378245, 298.3,         // Эллипсоид Красовского
                             6378137, 298.257223563, // Эллипсоид ПЗ-90.02 идентичен WGS-84 
                             direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378137, 298.257223563, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        /// <summary>
        /// Перевод СК-95 в ПЗ-90
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_268", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromCK95ToPZ90(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                             25.9, -130.94, -81.76,
                             0, 0, 0, 0,
                             6378245, 298.3,         //СК-95 использует эллипсоид Красовского
                             6378136, 298.257839303, //Эллипсоид ПЗ-90
                             direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378136, 298.257839303, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }


        /// <summary>
        /// Перевод СК-95 в ПЗ-90.02
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_269", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromCK95ToPZ9002(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                             24.83, -130.97, -81.74,
                             0, 0, -0.13, -0.22E-06,
                             6378245, 298.3,         // СК-95 использует эллипсоид Красовского
                             6378137, 298.257223563, // Эллипсоид ПЗ-90.02 идентичен WGS-84 
                             direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378137, 298.257223563, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        /// <summary>
        /// Перевод ПЗ-90 в WGS-84
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_270", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromPZ90ToWGS84(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                             -1.10, -0.30, -0.9,
                             0, 0, -0.20, -0.12E-06,
                             6378136, 298.257839303, // Эллипсоид ПЗ-90
                             6378137, 298.257223563, // Эллипсоид WGS-84
                             direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378137, 298.257223563, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        /// <summary>
        /// Перевод ПЗ-90.02 в WGS-84
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_271", typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromPZ9002ToWGS84(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                           -0.36, 0.08, 0.18,
                           0, 0, 0, 0,
                           6378137, 298.257223563, // Эллипсоид ПЗ-90.02 идентичен WGS-84
                           6378137, 298.257223563, // Эллипсоид WGS-84 
                           direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378137, 298.257223563, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        /// <summary>
        /// Перевод ПЗ-90 в ПЗ-90.02
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        [ConverterLocalized("Const_272",typeof(TimeZYX.Model.Properties.Resources))]
        public void ConvertFromPZ90ToPZ9002(ICollection<WellHeadGeoInfo> wellHeadGeoInfo, bool direct)
        {
            foreach (WellHeadGeoInfo geoInfo in wellHeadGeoInfo)
            {
                double w = geoInfo.WidthDec;
                double l = geoInfo.LengthDec;
                double h = geoInfo.Altitude;
                ConvertCoord(ref w, ref l, ref h,
                        1.07, 0.03, -0.02,
                        0, 0, 0.13, 0.22E-06,
                        6378136, 298.257839303, //Эллипсоид ПЗ-90
                        6378137, 298.257223563, //Эллипсоид ПЗ-90.02 идентичен WGS-84
                        direct);
                double x, y, z;
                BLHtoXYZ(w, l, h, 6378137, 298.257223563, out x, out y, out z);
                geoInfo.WidthDec = w; geoInfo.LengthDec = l; geoInfo.Altitude = h;
                geoInfo.X = x; geoInfo.Y = y; geoInfo.Z = z;
            }
        }

        #endregion

        #region Private methods


        /// <summary>
        /// Основной метод конвертации
        /// </summary>
        /// <param name="lat">Широта, десятичные градусы</param>
        /// <param name="lng">Долгота, десятичные градусы</param>
        /// <param name="h">Высота, м</param>
        /// <param name="dX">DX, линейные элементы трансформирования, в метрах</param>
        /// <param name="dY">DY, линейные элементы трансформирования, в метрах</param>
        /// <param name="dZ">DZ, линейные элементы трансформирования, в метрах</param>
        /// <param name="wX">WX, угловые элементы трансформирования, в секундах</param>
        /// <param name="wY">WY, угловые элементы трансформирования, в секундах</param>
        /// <param name="wZ">WZ, угловые элементы трансформирования, в секундах</param>
        /// <param name="mS">Дифференциальное различие масштабов</param>
        /// <param name="a1">Большая полуось первого эллипсоида, м</param>
        /// <param name="f1">Сжатие первого эллипсоида</param>
        /// <param name="a2">Большая полуось второго эллипсоида, м</param>
        /// <param name="f2">Сжатие второго эллипсоида</param>
        /// <param name="pr">Преобразование true - прямое, false - обратное</param>
        private void ConvertCoord(ref double lat, ref double lng, ref double h,
                                 double dX, double dY, double dZ,
                                 double wX, double wY, double wZ, double mS,
                                 double a1, double f1,
                                 double a2, double f2,
                                 bool pr)
        {

            aP = a1;
            alP = 1 / f1;

            aW = a2;
            alW = 1 / f2;

            dx = dX; dy = dY; dz = dZ;
            wx = wX; wy = wY; wz = wZ; 
            ms = mS;

            B = DegreesToRad(lat);
            L = DegreesToRad(lng);

            Init();

            int zn = 1;
            if (!pr) zn = -1;
            lat = lat + zn * dB(lat, lng, h) / 3600;
            lng = lng + zn * dL(lat, lng, h) / 3600;
            h = h + zn * dAlt(lat, lng, h);
        }


        private static double DegreesToRad(double val)
        {
            return val * Math.PI / 180;
        }

        private void Init()
        {
            e2P = 2 * alP - alP * alP;
            e2W = 2 * alW - alW * alW;

            a = (aP + aW) / 2;
            e2 = (e2P + e2W) / 2;
            da = aW - aP;
            de2 = e2W - e2P;

            Sin2B = Math.Pow(Math.Sin(B), 2);
            N = a * Math.Pow((1 - e2 * Sin2B), -0.5);
            M = a * (1 - e2) / Math.Pow((1 - e2 * Sin2B), 1.5);
        }

        private double dB(double Bd, double Ld, double H)
        {
            double _dB = ro/(M + H)*(N/a*e2*Math.Sin(B)*Math.Cos(B)*da +
                                     ((N*N)/(a*a) + 1)*N*Math.Sin(B)*Math.Cos(B)*de2/2 -
                                     (dx*Math.Cos(L) + dy*Math.Sin(L))*Math.Sin(B) + dz*Math.Cos(B)) -
                         wx*Math.Sin(L)*(1 + e2*Math.Cos(2*B)) +
                         wy*Math.Cos(L)*(1 + e2*Math.Cos(2*B)) -
                         ro*ms*e2*Math.Sin(B)*Math.Cos(B);
            return _dB;
        }

        private double dL(double Bd, double Ld, double H)
        {
            double _dL = ro/((N + H)*Math.Cos(B))*(-dx*Math.Sin(L) + dy*Math.Cos(L)) +
                         Math.Tan(B)*(1 - e2)*(wx*Math.Cos(L) + wy*Math.Sin(L)) - wz;
            return _dL;
        }

        private double dAlt(double Bd, double Ld, double H)
        {
            double dH = -a/N*da + N*Sin2B*de2/2 +
                        (dx*Math.Cos(L) + dy*Math.Sin(L))*Math.Cos(B) + dz*Math.Sin(B) -
                        N*e2*Math.Sin(B)*Math.Cos(B)*(wx/ro*Math.Sin(L) - wy/ro*Math.Cos(L)) +
                        ((a*a)/N + H)*ms;
            return dH;
        }

        private void BLHtoXYZ(double lat, double lng, double H, 
                              double a1, double f1,
                              out double x, out double y, out double z)
        {
            a = a1;
            alP = 1 / f1;

            B = DegreesToRad(lat);
            L = DegreesToRad(lng);

            Sin2B = Math.Pow(Math.Sin(B), 2);
            e2 = 2 * alP - alP * alP;
            N = a * Math.Pow((1 - e2 * Sin2B), -0.5);

            x = (N + H) * Math.Cos(B) * Math.Cos(L);
            y = (N + H) * Math.Cos(B) * Math.Sin(L);
            z = ((1 - e2)* N + H)* Math.Sin(B);

        }

        #endregion

        #region Fields

        const double ro = 206264.8062;

        // Первый эллипсоид 
        double aP, alP, e2P;

        // Второй эллипсоид 
        double aW, alW, e2W;

        // элементы трансформирования
        double dx, dy, dz, wx, wy, wz, ms;

        //Вспомогательные значения для преобразования эллипсоидов
        double a, e2, da, de2;

        //Вспомогательные значения
        double B, L, M, N, Sin2B;

        #endregion
    }

    [AttributeUsage(AttributeTargets.Method, AllowMultiple = false, Inherited = false)]
    public class ConverterAttribute : Attribute
    {
        private string name;
        public virtual string Name { get { return name; } }
        public MethodInfo MethodInfo { get; set; }

        public ConverterAttribute(string name)
        {
            this.name = name;
        }

        public override string ToString()
        {
            return name;
        }
    }

    [AttributeUsage(AttributeTargets.Method, AllowMultiple = false, Inherited = false)]
    public class ConverterLocalizedAttribute : ConverterAttribute
    {
        private readonly string resName;
        private readonly ResourceManager rm;

        /// <summary>
        /// Конструктор
        /// </summary>
        /// <param name="name">Наименование строки ресурса</param>
        public ConverterLocalizedAttribute(string name)
            : base(name)
        {
            resName = name;
            rm = Resources.ResourceManager;
        }

        public ConverterLocalizedAttribute(string name, Type type)
            : base(name)
        {
            resName = name;
            rm = new ResourceManager(type);
        }

        public override string Name
        {
            get
            {
                string str = rm.GetString(resName);
                return (!string.IsNullOrEmpty(str)) ? str : resName;
            }
        }
    }
}
