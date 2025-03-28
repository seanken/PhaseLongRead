package AlleleMASSeq;
import java.io.FileWriter;
import java.lang.*;
import java.io.*;
import htsjdk.samtools.*;
import java.util.*;
import java.util.Map.Entry;
import java.util.Map;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;


public class Long_ASE_Count
{

    protected SAMFileWriter savFil1, savFil2, savFil3;
    protected File bamFile, outBam1, outBam2, outBam3;
    protected HashMap<SNPVal, Boolean> map;
    protected HashMap<Integer, Boolean> posMap;
    //protected HashMap<String, String> UMIMap;
    protected HashMap<String,Integer> counts;
    protected int Good;
    protected int Bad;
    protected int numRem;
    //protected String gene_tag;
    //protected String umi_tag;

    public Long_ASE_Count(String args[]) throws Exception
    {
        //ArrayList<String[]> ret=new ArrayList<String[]>();
        this.bamFile = new File(args[0]); 
        File vcfFile=new File(args[1]);
        this.outBam1 = new File(args[2]);
        this.outBam2 = new File(args[3]);
        this.outBam3 = new File(args[4]);
        //int ind=Integer.parseInt(args[2]);
        int ind=0; //assums vcf is prepared 
        String saveFile=args[2];
        //this.gene_tag=args[3];
        //this.umi_tag=args[4];


        

        print("Create VCF HashMap");

        //this.map= 
        HashMapVCF(vcfFile,ind);
        System.out.println(this.map.size());
    }


    public void IterateOverAll() throws Exception
    {

        print("Read in long read data!");
        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
        SAMRecordIterator r = sr.iterator();
        this.savFil1=new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),true, this.outBam1);
        this.savFil2=new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),true, this.outBam2);
        this.savFil3=new SAMFileWriterFactory().makeSAMOrBAMWriter(sr.getFileHeader(),true, this.outBam3);
        int iter=0;
        this.Good=0;
        this.Bad=0;
        this.numRem=0;

        //this.UMIMap=new HashMap<String, String>(); 

        while(r.hasNext()) {
            iter=iter+1;

            if(iter%10000==0)
            {
                System.out.println(this.Good);
                System.out.println(this.Bad);
                System.out.println(this.numRem);
                System.out.println((float)this.Good/(float)(this.Good+this.Bad));
                System.out.println((float)this.Good/(float)(iter));
                System.out.println(iter);
            }

            try{
                SAMRecord read=r.next();
                this.ProcessRead(read,iter);
            }catch(Exception e){
                print("Yuck");
                break;
            }

        }
        
        r.close();
        sr.close();
        this.savFil1.close();
        this.savFil2.close();
        this.savFil3.close();

    }




    public void HashMapVCF(File vcfFile,int ind) throws Exception
    {

     
        Scanner sc = new Scanner(vcfFile); 

        HashMap<SNPVal,Boolean> vcfMap=new HashMap<SNPVal,Boolean>();
        this.posMap=new HashMap<Integer,Boolean>();

        int iter=0;

        while (sc.hasNextLine()) 
        {
            iter=iter+1;

            if(iter %100000==0)
            {
                System.out.println(iter);
            }

            String line=sc.nextLine();
            if('#'==line.charAt(0))
            {
                continue;
            }

            String[] split_line=line.split("\\s+");


            String chrom=split_line[0];
            String pos=split_line[1];
            String ref=split_line[3];
            String alt=split_line[4];
            int genoPos=9+ind;
            String geno=split_line[genoPos];
            String[] split_geno=geno.split(":");
            geno=split_geno[0];

            if(ref.length()>1 || alt.length()>1)
            {
                continue;
            }
                
            int posInt=Integer.parseInt(pos);
            chrom=chrom.replace("chr","");

            SNPVal key1=new SNPVal(posInt,chrom,ref.charAt(0));
            SNPVal key2=new SNPVal(posInt,chrom,alt.charAt(0));

            Boolean val1=true;
            Boolean val2=false;

            //int posInt=Integer.parseInt(pos);

            if(geno.equals("0|1") || geno.equals("0/1"))
            {
                vcfMap.put(key1,val1);
                vcfMap.put(key2,val2);
                this.posMap.put(posInt,val1);
            }


            if(iter %100000==0)
            {
                System.out.println(key1);
                System.out.println("");
            }




            if(geno.equals("1|0") || geno.equals("1/0"))
            {
                vcfMap.put(key1,val2);
                vcfMap.put(key2,val1);
                this.posMap.put(posInt,val1);
            }


        }

        this.map=vcfMap;

        //return(vcfMap);

    }



    public void ProcessRead(SAMRecord read,int inter)
    {
        String seq=read.getReadString();

        //Character isprim=read.getCharacterAttribute("tp");
        String suppmaps=read.getStringAttribute("SA");

        if(!(suppmaps==null))
        {
            return;
        }

        //Add in multithis.map filtering

        
        this.numRem=this.numRem+1;

        //getReferencePositionAtReadPosition
        int start=read.getStart();
        String chrom=read.getContig();
        chrom=chrom.replace("chr","");
        int readLen=read.getReadLength();
        int All1=0;
        int All2=0;

        int[] posLst=getPosition(read,readLen,start);

        for(int i=0;i<readLen;i=i+1)
        {

            //int pos=read.getReferencePositionAtReadPosition(i+1);
            int pos=posLst[i];
            //pos=pos;
            //Just to double check our code
            if(inter<20)
            {
                int pos2=read.getReferencePositionAtReadPosition(i+1);
                if(pos!=pos2)
                {
                    print("Error with positon!");
                }
            }

            if(pos<1){continue;}

            if(!this.posMap.containsKey(pos))
            {
                continue;
            }


            //String base=Character.toString(seq.charAt(i));
            //String base=Character.toString(seq.charAt(i));
            char base=seq.charAt(i);

            //String key=chrom+"_"+Integer.toString(pos)+"_"+base;
            SNPVal key=new SNPVal(pos,chrom,base);


            if(this.map.containsKey(key))
            {
                Boolean val=this.map.get(key);
                if(val)
                {
                    All1=All1+1;
                }
                else
                {
                    All2=All2+1;
                }

            }


        }

        int Tot=All1+All2;
        if(Tot<1){return;}
        float rat=(float)Math.max(All1,All2)/(float)Tot;
        read.setAttribute("AL",All1);
        read.setAttribute("AM",All2);

        if(rat>.95)
        {
            this.Good=this.Good+1;
        }
        else
        {
            this.Bad=this.Bad+1;
            this.savFil3.addAlignment(read);
            return;
        }



        String val="None";
        if(All1>All2)
        {
            val="All1";
            this.savFil1.addAlignment(read);
        }
        if(All2>All1)
        {
            val="All2";
            this.savFil2.addAlignment(read);
        }

        read.setAttribute("AL",val);




    }



    public static void print(String toPrint)
    {

    System.out.println(toPrint);

    }



    //Gets an array with position information (use this to speed up process)
    public static int[] getPosition(SAMRecord read,int readLen,int start) //throws Exception
    {
        int curPosRead=0;
        int curPosRef=start;
        int[] posLst=new int[readLen];
        List<CigarElement> cigar=read.getCigar().getCigarElements();
        for (int i = 0; i < cigar.size(); i++) {
            int len=cigar.get(i).getLength();
            CigarOperator op=cigar.get(i).getOperator();
            if(op.consumesReadBases() && op.consumesReferenceBases())
            {
                for(int j=0;j<len;j++)
                {
                    posLst[curPosRead]=curPosRef;
                    curPosRead=curPosRead+1;
                    curPosRef=curPosRef+1;
                }
            }
            else if(op.consumesReadBases())
            {
                for(int j=0;j<len;j++)
                {
                    posLst[curPosRead]=0;
                    curPosRead=curPosRead+1;
                }
            }
            else if(op.consumesReferenceBases())
            {
                curPosRef=curPosRef+len;
            }
        }

        return(posLst);
        
    }


    public static HashMap<String,Integer> HashMapVCF_SV(File vcfFile,int ind) throws Exception
    {
     
        Scanner sc = new Scanner(vcfFile); 

        HashMap<String,Integer> vcfMap=new HashMap<String,Integer>();

        int iter=0;

        while (sc.hasNextLine()) 
        {
            iter=iter+1;

            if(iter %100000==0)
            {
                System.out.println(iter);
            }

            String line=sc.nextLine();
            if('#'==line.charAt(0))
            {
                continue;
            }

            String[] split_line=line.split("\\s+");


            String chrom=split_line[0];
            String pos=split_line[1];
            String ref=split_line[3];
            String alt=split_line[4];
            int genoPos=9+ind;
            String geno=split_line[genoPos];
            String[] split_geno=geno.split(":");
            geno=split_geno[0];

            

            chrom=chrom.replace("chr","");

            String key=chrom+"_"+pos;
            int val=1;

            if(geno.equals("0|1") || geno.equals("0/1"))
            {
                vcfMap.put(key,val);
            }


            if(iter %100000==0)
            {
                System.out.println(key);
                System.out.println("");
            }




            if(geno.equals("1|0") || geno.equals("1/0"))
            {
                vcfMap.put(key,val);
            }


        }

        return(vcfMap);

    }
}

