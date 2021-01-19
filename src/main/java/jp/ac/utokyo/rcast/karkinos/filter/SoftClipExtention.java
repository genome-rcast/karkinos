/*
Copyright Hiroki Ueda

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
package jp.ac.utokyo.rcast.karkinos.filter;

import jp.ac.utokyo.rcast.karkinos.utils.TwoBitGenomeReader;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;

public class SoftClipExtention {

    //add 2020/12/17 for FFPE anneling near repeat, extends soft clipping
    public static void extendSoftclip(SAMRecord sam, TwoBitGenomeReader tgr) {
        try {
            if(sam.getCigar().getCigarElements().size()>=2) {
                leftExtention(sam, tgr);
                rightExtention(sam, tgr);
            }
        }catch(Exception ex){
            //could not extend soft clips
            ex.printStackTrace();
        }
    }

    private static void leftExtention(SAMRecord sam, TwoBitGenomeReader tgr) {
        List<CigarElement> list = sam.getCigar().getCigarElements();
        CigarElement first = list.get(0);
        CigarElement second = list.get(1);
        int extraCheckLen = 20;

        if(isClipped(first)&&isM(second)&&second.getLength()>extraCheckLen){

            int alignmentStart = sam.getAlignmentStart();
            int miscount = 0;
            int lastidx = 0;
            int m = 0;
            int clipidx = 0;
            for(int n=alignmentStart;n<=alignmentStart+extraCheckLen;n++){
                //
                char refNuc = tgr._getGenomeNuc(n,true);
                int idx = sam.getReadPositionAtReferencePosition(n)-1;
                if(idx>=0 && idx<sam.getReadLength()){
                    char readNuc = sam.getReadString().charAt(idx);
                    if(!equalnuc(refNuc,readNuc)){
                        miscount++;
                        lastidx = idx;
                    }
                }
                if(m==9){
                    if(miscount>0){
                        clipidx = lastidx;
                    }
                }else if(m==14){
                    if(miscount>1){
                        clipidx = lastidx;
                    }
                }else if(m==19){
                    if(miscount>2){
                        clipidx = lastidx;
                    }
                }
                m++;
            }
            //
            if(clipidx>0){
                Cigar cgn = new Cigar();
                CigarElement sc = new CigarElement(clipidx+1,CigarOperator.S);
                cgn.add(sc);
                CigarElement mc = new CigarElement((first.getLength()+second.getLength())-(clipidx+1),CigarOperator.M);
                cgn.add(mc);
                for(m=2;m<list.size();m++){
                    cgn.add(list.get(m));
                }
                sam.setCigar(cgn);
                int alstart = sam.getStart()+(clipidx+1-first.getLength());
                sam.setAlignmentStart(alstart);
                int  nm = sam.getIntegerAttribute("NM");
                nm = nm - miscount;
                sam.setAttribute("NM",nm);
            }
        }
    }

    private static void rightExtention(SAMRecord sam, TwoBitGenomeReader tgr) {
        List<CigarElement> list = sam.getCigar().getCigarElements();
        int size = list.size();
        CigarElement first = list.get(size-1);
        CigarElement second = list.get(size-2);
        int extraCheckLen = 20;

        if(isClipped(first)&&isM(second)&&second.getLength()>extraCheckLen){

            int alignmentEnd = sam.getAlignmentEnd();
            int miscount = 0;
            int lastidx = 0;
            int m = 0;
            int clipidx = 0;
            for(int n=alignmentEnd;n>=alignmentEnd-extraCheckLen;n--){
                //
                char refNuc = tgr._getGenomeNuc(n,true);
                int idx = sam.getReadPositionAtReferencePosition(n)-1;
                if(idx>=0 && idx<sam.getReadLength()){
                    char readNuc = sam.getReadString().charAt(idx);
                    if(!equalnuc(refNuc,readNuc)){
                        miscount++;
                        lastidx = idx;
                    }
                }
                if(m==9){
                    if(miscount>0){
                        clipidx = lastidx;
                    }
                }else if(m==14){
                    if(miscount>1){
                        clipidx = lastidx;
                    }
                }else if(m==19){
                    if(miscount>2){
                        clipidx = lastidx;
                    }
                }
                m++;
            }
            //
            if(clipidx>0){
                //
                Cigar cgn = new Cigar();
                for(m=0;m<list.size()-2;m++){
                    cgn.add(list.get(m));
                }
                int sclen = sam.getReadLength()-(clipidx+1);
                int mlen = second.getLength() - (sclen - first.getLength());
                CigarElement mc = new CigarElement(mlen,CigarOperator.M);
                CigarElement sc = new CigarElement(sclen,CigarOperator.S);
                cgn.add(mc);
                cgn.add(sc);
                sam.setCigar(cgn);
                int  nm = sam.getIntegerAttribute("NM");
                nm = nm - miscount;
                sam.setAttribute("NM",nm);
            }
        }
    }

    private static boolean equalnuc(char refNuc,char readNuc) {
        return Character.toUpperCase(refNuc) == Character.toUpperCase(readNuc);
    }

    private static boolean isClipped(CigarElement ce) {
        return ce.getOperator().equals(CigarOperator.S);
    }

    private static boolean isM(CigarElement ce) {
        return ce.getOperator().equals(CigarOperator.M);
    }
}
