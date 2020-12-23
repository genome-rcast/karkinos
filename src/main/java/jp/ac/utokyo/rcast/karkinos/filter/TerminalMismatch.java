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

import java.util.List;

public class TerminalMismatch {
    //add 2020/12/17 for FFPE anneling near repeat, extends soft clipping
    public static int terminalMismatch(SAMRecord sam, TwoBitGenomeReader tgr, int extraCheckLen) {
        try {
            //does not check Indel read
            if(containIndel(sam)){
                return 0;
            }
            return Math.max(leftMismatch(sam,tgr,extraCheckLen),rightMismatch(sam,tgr,extraCheckLen));

        }catch(Exception ex){
            //
            ex.printStackTrace();
        }
        return 0;
    }

    private static boolean containIndel(SAMRecord sam){
        List<CigarElement> l = sam.getCigar().getCigarElements();
        for(CigarElement ce: l){
            if(ce.getOperator().equals(CigarOperator.I))return true;
            if(ce.getOperator().equals(CigarOperator.D))return true;
        }
        return false;
    }

    private static int leftMismatch(SAMRecord sam, TwoBitGenomeReader tgr, int extraCheckLen) {
        int alignmentStart = sam.getAlignmentStart();
        int miscount = 0;
        CigarElement first = sam.getCigar().getCigarElements().get(0);
        int scidx = 0;
        if(isClipped(first)){
            scidx = first.getLength();
        }
        for(int n=alignmentStart;n<=alignmentStart+extraCheckLen;n++){
            //
            char refNuc = tgr._getGenomeNuc(n,true);
            int idx = sam.getReadPositionAtReferencePosition(n)-1+scidx;
            if(idx>=0 && idx<sam.getReadLength()){
                char readNuc = sam.getReadString().charAt(idx);
                if(!equalnuc(refNuc,readNuc)){
                    miscount++;
                }
            }
        }
        return miscount;
    }

    private static int rightMismatch(SAMRecord sam, TwoBitGenomeReader tgr, int extraCheckLen) {
        int alignmentEnd = sam.getAlignmentEnd();
        int miscount = 0;
        CigarElement first = sam.getCigar().getCigarElements().get(0);
        int scidx = 0;
        if(isClipped(first)){
            scidx = first.getLength();
        }
        for(int n=alignmentEnd;n>=alignmentEnd-extraCheckLen;n--){
            //
            char refNuc = tgr._getGenomeNuc(n,true);
            int idx = sam.getReadPositionAtReferencePosition(n)-1+scidx;
            if(idx>=0 && idx<sam.getReadLength()){
                char readNuc = sam.getReadString().charAt(idx);
                if(!equalnuc(refNuc,readNuc)){
                    miscount++;
                }
            }
        }
        return miscount;
    }

    private static boolean equalnuc(char refNuc,char readNuc) {
        return Character.toUpperCase(refNuc) == Character.toUpperCase(readNuc);
    }

    private static boolean isClipped(CigarElement ce) {
        return ce.getOperator().equals(CigarOperator.S);
    }
}
