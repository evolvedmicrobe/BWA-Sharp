using System;

namespace ccs_caller
{
    public class ZMWInfo
    {
        public string Movie;
        public int ZMW;
        public float PredictedRawAccuracy;
        public float SnrC;
        public float SnrG;
        public float SnrT;
        public float SnrA;
        public ushort Length;
        public ushort NumPasses;
        public float PredictedCCSAccuracy;
        public ZMWInfo (string line)
        {
            var sp = line.Trim().Split (',');
            Movie = String.Intern (sp [0]);
            ZMW = Convert.ToInt32 (sp [1]);
            PredictedRawAccuracy = Convert.ToSingle (sp [2]);
            SnrT = Convert.ToSingle (sp [3]);
            SnrG = Convert.ToSingle (sp [4]);
            SnrA = Convert.ToSingle (sp [5]);
            SnrC = Convert.ToSingle (sp [6]);
            Length = Convert.ToUInt16 (sp [7]);
            NumPasses = Convert.ToUInt16 (sp [8]);
            PredictedCCSAccuracy = Convert.ToSingle (sp [9]);
        }
    }
}

