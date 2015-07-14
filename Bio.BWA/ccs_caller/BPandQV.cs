using System;

namespace ccs_caller
{
    public struct BPandQV
    {
        public readonly byte BP;
        public readonly byte QV;
        public BPandQV (byte bp, byte qv)
        {
            BP = bp;
            QV = qv;
        }
    }
}

