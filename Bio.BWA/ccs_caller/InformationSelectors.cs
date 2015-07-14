using System;
using System.Collections.Generic;
using System.Linq;
using Bio;
namespace ccs_caller
{
    public static class InformationSelectors
    {
        public class InfoSelector {
            public Func<AlignedCCSRead, Variant, object> Function;
            public string Name;
            public InfoSelector(string name, Func<AlignedCCSRead, Variant, object> f) {
                Function = f;
                Name = name;
            }
        }

        public static string GetHeader()
        {
            return String.Join("," , DataGetters.Select(x => x.Name).ToArray());
        }

        public static string ConvertObject(object o) {
            if (o == null) {
                return "NA";
            } else {
                return o.ToString ();
            }
        }

        public static string GetDataLine(Variant v, AlignedCCSRead read) {
            return String.Join("," , DataGetters.Select(x => ConvertObject(x.Function(read, v))).ToArray());
        }


        public static List<InfoSelector> DataGetters = new List<InfoSelector> () {
            new InfoSelector ("ZMW", (x, y) => x.Zmw.ZMW),
            new InfoSelector ("Movie", (x, y) => x.Zmw.Movie),
            new InfoSelector ("NumPasses", (x, y) => x.Zmw.NumPasses),
            new InfoSelector ("RevComp", (x, y) => x.OriginallyReverseComplemented),
            new InfoSelector ("MapQV", (x, y) => x.MapQV),
            new InfoSelector ("AlignedLength", (x, y) => x.End - x.Start),
            new InfoSelector ("SnrC", (x, y) => x.Zmw.SnrC),
            new InfoSelector ("SnrG", (x, y) => x.Zmw.SnrG),    
            new InfoSelector ("SnrT", (x, y) => x.Zmw.SnrT),
            new InfoSelector ("SnrA", (x, y) => x.Zmw.SnrA),
            new InfoSelector ("Pos", (x,y) => y.StartPosition),
            new InfoSelector ("Ref", (x,y) => y.RefName),
            new InfoSelector ("Type", (x,y) => y.Type),
            new InfoSelector("NumVariantsInRead", (x,y) => x.Variants.Count),
            new InfoSelector ("IsVariant,QV,Seq", GetVariantInformation),
            new InfoSelector("InHP,HPLen,HPBase,DelorIns", (x,y) => { var q = (y as IndelVariant);
                if(q == null) 
                    return "NA,NA,NA, NA"; 
                else 
                    return q.InHomopolymer + "," + q.HomopolymerLengthInReference + "," + q.HomopolymerBase + "," + q.InsertionOrDeletion.ToString();}) 

        };

        public static string GetVariantLine(AlignedCCSRead read, Variant v)
        {
            var arr = DataGetters.Select (x => ConvertObject (x.Function (read, v))).ToArray();
            return String.Join (",", arr );
        }

        // TODO: Method needs verification and work, the relevant QV is not well defined for indels and I 
        // don't want to handle the complex (and typically totally irrelevant) full logic here.
        private static string GetVariantInformation(AlignedCCSRead read, Variant variant) {
            var pos = variant.StartPosition;
            var current = read.Start;
            var refSeq = read.AlignedReferenceSeq;
            int variant_position_in_alignment = -1;
            for (int i = 0; i < refSeq.Length; i++) {
                if (refSeq [i] != (byte)'-') {
                    if (current == pos) {
                        variant_position_in_alignment = i;
                        break;
                    }
                    current++;
                }
            }

            // Now what type of variant are we dealing with?
            bool isVariant = true;
            byte QV = 255;
            string seq = String.Empty;
            var cq = read.AlignedQuerySeq;
            var cr = read.AlignedReferenceSeq;
            if (variant.Type == VariantType.SNP) {
                var cur = cq[variant_position_in_alignment];
                isVariant = refSeq [variant_position_in_alignment] != cur.BP;
                QV = cur.QV;
                seq = Convert.ToString ((char)cur.BP);
            } else if (variant.Type == VariantType.INDEL) {
                var nv = variant as IndelVariant;

                if (nv.InsertionOrDeletion == IndelType.Insertion) {
                    // This implies that the next reference position should be a ref.
                    // TODO: Account for the fact that the next position might not be there
                    isVariant = cr [variant_position_in_alignment + 1] == DnaAlphabet.Instance.Gap;
                    // Now to get the QV value 
                    if (isVariant) {
                        List<char> chars = new List<char> ();
                        int cp = variant_position_in_alignment + 1;
                        while (cr [cp] == DnaAlphabet.Instance.Gap) {
                            chars.Add ((char)cq [cp].BP);
                            cp++;
                        }
                        seq = new string (chars.ToArray ());
                    } else {
                        seq = "NA";
                    }
                    if (!read.OriginallyReverseComplemented) {
                        // This should represent the probability of an insertion,
                        // or a deletion prior to the current base, so in either case should be the QV for this
                        QV = cq [variant_position_in_alignment + 1].QV;
                    } else {
                        // We need to account for any homopolymers, as all the deletions will be on that
                        int cp = variant_position_in_alignment + 1;
                        var qbp = cq [cp].BP;
                        // If at the end, we have to assume the QV is there
                        // this is an edge condition that is bad though.
                        while ((cp+1) != cq.Length && cq [cp + 1].BP == qbp) {
                            cp++;
                        }

                        QV = cq [cp].QV;
                    }
                } else if (nv.InsertionOrDeletion == IndelType.Deletion) {
                    isVariant = cq [variant_position_in_alignment + 1].BP == DnaAlphabet.Instance.Gap;
                    // Now to get the QV value 
                    int cp = variant_position_in_alignment + 1;
                    if (isVariant) {
                        List<char> chars = new List<char> ();
                        while (cq [cp].BP == DnaAlphabet.Instance.Gap) {
                            chars.Add ((char)cr [cp]);
                            cp++;
                        }
                        seq = new string (chars.ToArray ());
                    } else {
                        seq = "NA";
                    }
                    if (!read.OriginallyReverseComplemented) {
                        // This should represent the probability of a deletion prior to the current base, 
                        QV = isVariant ? cq [cp].QV : cq [variant_position_in_alignment].QV;
                    } else {
                        // We need to account for any homopolymers, as all the deletions will be on that
                        var qbp = cq [cp].BP;
                        while ((cp+1)  < cq.Length && cq [cp + 1].BP == qbp) {
                            cp++;
                        }
                        QV = cq [cp].QV;
                    }
                }

            }
            return isVariant + "," + QV + "," + seq;
        }
    }
}

