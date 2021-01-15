public class TEST_RNA_Hybridization {
  
  public static void main(String[] args) {
    
    if (args.length < 2) {
      System.err.println("\nTo execute this program, two command line arguments are required corresponding to the names of two FASTA files, each containing an RNA sequence. For example,\n");
      System.err.println("\tjava TEST_RNA_Hybridization sequences/RNA_2a.txt sequences/RNA_2b.txt\n");
      return;
    }

    RNA_Hybridization_Sols r = new RNA_Hybridization_Sols(args[0], args[1]);
    r.computeHybridizationWithMinimumEnergy();
    r.determineHybridization();
    System.out.println(r.arrayToString());
    System.out.println("\nHybridization with lowest energy:\n");
    System.out.println("\t" + r.get_s1_hybridization());
    System.out.println("\t" + r.get_hybridization());
    System.out.println("\t" + r.get_s2_hybridization());
    System.out.println("\nOptimal hybridization energy: " + r.get_optimalEnergy() + "\n");    
  }  
  
}

