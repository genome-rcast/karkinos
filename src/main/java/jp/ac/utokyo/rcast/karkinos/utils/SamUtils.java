package jp.ac.utokyo.rcast.karkinos.utils;


import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SamUtils {

    public static boolean lowmap(SAMRecord sam) {

        if(sam.getReadUnmappedFlag()
         || (sam.getMappingQuality() <= 3)
         || sam.getNotPrimaryAlignmentFlag()){
            return true;
        }


        int seqlen = sam.getReadLength();
        int match = 0;
        int indelcnt = 0;
        int indellen = 0;
        int sclen = 0;
        for (CigarElement ce : sam.getCigar().getCigarElements()) {
            // if(ce.getOperator() == CigarOperator.S)return true;
            if (ce.getOperator() == CigarOperator.I) {
                indelcnt++;
                indellen = indellen + ce.getLength();
            }
            if (ce.getOperator() == CigarOperator.D) {
                indelcnt++;
                indellen = indellen + ce.getLength();
            }
            if (ce.getOperator() == CigarOperator.S) {
                sclen = sclen + ce.getLength();
            }

        }

        if (indelcnt >= 3)return true;

        Integer nm = sam.getIntegerAttribute("NM");
        if(nm!=null){
            nm = nm - indellen;
        }else{
            nm=0;
        }
        if (sam.getMappingQuality() <= 15 && nm>=2) {
            return true;
        }
        double errorratio = (double)nm/(double)(seqlen-sclen);
        if(errorratio>0.15){
            return true;
        }

        return false;
    }



}
