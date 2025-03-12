package AlleleMASSeq;

public class SNPVal
{
    protected int pos;
    protected String chrom;
    protected char base;

    public SNPVal(int pos, String chrom, char base) 
    {
        this.pos=pos;
        this.base=base;
        this.chrom=chrom;
    }

    @Override public boolean equals(Object o)
    {
        if (this.getClass() != o.getClass()) {
            return false;
        }
        SNPVal snp=(SNPVal) o;
        boolean isSame=(this.pos==snp.getPos() && this.base==snp.getBase() && this.chrom.equals(snp.getChrom()));
        return(isSame);
    }

    public int getPos()
    {
        return(this.pos);
    }

    public char getBase()
    {
        return(this.base);
    }

    public String getChrom()
    {
        return(this.chrom);
    }

    @Override public int hashCode() 
    {
        return java.util.Objects.hash(this.pos,this.chrom,this.base);
    }

    @Override public String toString() 
    {
        String val=this.chrom+"_"+String.valueOf(this.pos)+"_"+Character.toString(this.base);
        return(val);
    }

}