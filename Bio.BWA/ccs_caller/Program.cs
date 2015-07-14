using System;
using Bio.IO.FastQ;
using Bio;
using System.Linq;
using System.Collections.Generic;
using System.IO;


namespace ccs_caller
{
    using st = Tuple<string, string>;
    class MainClass
    {
        public static void Main (string[] args)
        {
            var fnames = new List<st> () {
                new st ("/Users/nigel/git/cafe-quality/NotTracked/CCS_Human/pulsewidth", "pulsewidth.csv"), 
                new st (    "/Users/nigel/git/cafe-quality/NotTracked/CCS_Human/master", "master.csv"),
                new st("/Users/nigel/git/cafe-quality/NotTracked/CCS_Human/mqv_deltag", "mqvdeltag.csv")
            };

            foreach(var t in fnames) {
                var direc = t.Item1;
            //var fname = "/Users/nigel/git/cafe-quality/NotTracked/CCS_Human/m150324_103244_42251_c100709512550000001823139404221516_s1_p0.1.ccs.fastq";

            // Load CSV Files 
            var zmwInfo = LoadZMWInfo( direc);

            var fastq = new FastQParser ();
                var fastqs = (new DirectoryInfo (direc)).GetFiles ().Where (z => z.Name.EndsWith ("1.ccs.fastq", StringComparison.Ordinal)).Select(p=>p.FullName);
                var seq = fastqs.SelectMany( fq =>  fastq.Parse (fq).Select (s =>  {
                    var r = s as QualitativeSequence;
                    var qid = r.ID;
                    var sp = qid.Split('/');
                    var movie = String.Intern(sp[0]);
                    var zmw = Convert.ToInt32(sp[1]);
                    var info = zmwInfo[movie][zmw];
                    var ccs=  CCSReadAligner.AlignRead(r, info.NumPasses);
                    if(ccs!=null) {
                        ccs.Zmw = info;
                    }
                    return ccs;
                })).ToList();
                CCSReadAligner.PrintRegionTree ("Tree.csv");    
            var totalReads = seq.Count;
            Console.WriteLine ("Total CCS Reads: " + totalReads);
            var seqf = seq.Where (o => o != null && o.Zmw.PredictedCCSAccuracy >= .999).ToList ();
            seq = null;
            var passFilters = seqf.Count;
            Console.WriteLine ("Passing Filters: " + passFilters + "  - " + ((passFilters / (double)totalReads).ToString ("P1")));
                var variants = seqf.SelectMany (p => p.Variants).Distinct ().Where(z => {
                    var reg = CCSReadAligner.GetRegionInformation(z);
                    return reg.ObservationCount > 10 && (reg.End-reg.Start) > 150;})
                    .ToList();
            OutputFile (variants, seqf, t.Item2);
            Console.WriteLine ("Total Variants: " + variants.Count);
            Console.WriteLine (seq);
            }
        }

        
        public static void OutputFile(List<Variant> variants, List<AlignedCCSRead> reads, string fname)
        {
            // Hack to avoid accounting for variants on the edges of alignments.
            const int EDGE_WINDOW = 5;
            variants.Sort ();
            reads.Sort ();

            StreamWriter sw = new StreamWriter(fname);
            sw.WriteLine(String.Join(",", "VariantNumber", InformationSelectors.GetHeader()));
            int currentVariant = 0;

            int firstReadInRegion = 0;
            foreach (var v in variants) 
            {
                // Let's get all the reads that overlap with this variant.
                List<AlignedCCSRead> cReads  = new List<AlignedCCSRead>();
                int current = firstReadInRegion;
                // This search is slow as hell, need to decide on the maximum window size so I can bound this.
                while (true) {
                    var cr = reads [current];
                    if ((cr.Start + EDGE_WINDOW) <= v.StartPosition && (cr.End - EDGE_WINDOW) > v.StartPosition) {
                        cReads.Add (cr);
                    } else if (String.CompareOrdinal(cr.Reference, v.RefName) > 0) {
                        // If we are on the next chromosome, this is the definite end
                        break;
                    } else if (cr.End < v.StartPosition && String.CompareOrdinal(v.RefName, cr.Reference) <= 0) {
                        firstReadInRegion = current;
                    }
                    current++;
                }

                // Now for each let's make a variant line 
                var refN = v.RefName;
                var pos = v.StartPosition.ToString();
                var lines = cReads.Select(x=> InformationSelectors.GetVariantLine(x, v));
                foreach(var l in lines)
                {
                    var line = l;
                    sw.WriteLine(currentVariant + ", " + line);
                }
                currentVariant++;
            }
            sw.Close ();

        }
        
        public static Dictionary<string, Dictionary<int, ZMWInfo>> LoadZMWInfo(string dirName) {
            DirectoryInfo di = new DirectoryInfo (dirName);
            var files = di.GetFiles ().Where (p => p.Name.EndsWith ("1.ccs.csv", StringComparison.Ordinal));
            var toReturn = new Dictionary<string, Dictionary<int, ZMWInfo>> ();
            foreach (var f in files) {
                var data = File.ReadAllLines (f.FullName).Skip (1).Select (k => new ZMWInfo (k));

                string movie = String.Empty;
                var dict = new Dictionary<int, ZMWInfo> ();
                bool movieSet = false;
                foreach (var d in data) {
                    if (movieSet) {
                        if (d.Movie != movie) {
                            throw new Exception ("There were two movies in the same ccs.csv file");
                        }
                    } else {
                        movie = d.Movie;
                        movieSet = true;
                    }
                    dict [d.ZMW] = d;
                }
                if (toReturn.ContainsKey (movie)) {
                    var odict = toReturn [movie];
                    foreach (var kv in dict) {
                        odict [kv.Key] = kv.Value;
                    }
                } else {
                    toReturn [movie] = dict;
                }
            }
            return toReturn;
        }

    }
}
