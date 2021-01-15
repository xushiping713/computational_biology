import java.util.*;
import java.io.*;

public class RNA_Hybridization {
  
  // Instance variables
  private String s1;
  private String s2;
  private String s2_reverse;
  private Energies e;
  private double[][] H;
  private double optimalEnergy;
  private int optimal_i;
  private int optimal_j;
  private String upper;
  private String hybridization;
  private String lower;
  
  // Constructor
  public RNA_Hybridization(String fileName1, String fileName2) {
    s1 = SequenceOps.sequenceFromFastaFile(fileName1).toUpperCase().replace('T', 'U');
    s2 = SequenceOps.sequenceFromFastaFile(fileName2).toUpperCase().replace('T', 'U');
    s2_reverse = new StringBuilder(s2).reverse().toString();
    e = new Energies("energyData", s1, s2_reverse);
    H = new double[s1.length()][s2.length()];
    optimalEnergy = 0.0;
    upper="";
    hybridization="";
    lower="";
  }
  
  // Public instance methods
  // For 2 RNA sequences, fill in the dynamic programming table, H, to compute the hybridization
  // that has the lowest energy
  public void computeHybridizationWithMinimumEnergy() {
    
    for (int i=0; i<s1.length(); i++) {
      for (int j=0; j<s2.length(); j++) {
      
        // Fill in H table entry
        if (j==0 || i==0) {
          H[i][j] = 0.0; // Fill in first column and row of H matrix with 0
        }
        else { // compute the minimum energy of hybridization between subsequence s1 and subsequence s2
          
          // stacking region
          double caseH_1 = e.stacking(i-1, j-1, i, j) + H[i-1][j-1];
          
          // bulge loop (first sequence)
          double caseH_2 = 0.0;  
          for (int x=0; x<i-1; x++) {
            double temp = e.bulge(x,j-1,i,j) + H[x][j-1];
            if (temp < caseH_2) caseH_2 = temp;
          }
          
          // bulge loop (second sequence)
          double caseH_3 = 0.0;
          for (int y=0; y<j-1; y++) {
            double temp = e.bulge(i-1,y,i,j)+H[i-1][y];
            if (temp < caseH_3) caseH_3 = temp;
          }
          // interior loop
          double caseH_4 = 0.0;
          for (int x=0; x<i-1; x++) {
            for (int y=0; y<j-1; y++) {
              double temp = e.interior(x,y,i,j) + H[x][y];
              if (temp <caseH_4) caseH_4 = temp;
            }
          }
          
          H[i][j] = min(caseH_1, caseH_2, caseH_3, caseH_4, 0);
          if (H[i][j]<optimalEnergy) { 
            optimal_i = i;
            optimal_j = j;
            optimalEnergy = H[i][j];
          }
        }
      }
    }
  }
  
  // returns energy of the optimal hybridization
  public double get_optimalEnergy() {
    return optimalEnergy; 
  }
  
  // Backtrack through dynamic programming tables to compute three lines that represents 
  // the hybirdization with optimal energy
  public void determineHybridization() {
    upper = s1.substring(optimal_i+1,s1.length());
    lower = s2_reverse.substring(optimal_j+1, s2.length());
    determineHybridization(optimal_i,optimal_j);
  }
  
  public void determineHybridization(int i, int j) {

    if (H[i][j]==0) {
      upper = spaces(Math.max(i,j)-i)+s1.substring(0,i+1)+upper;
      hybridization = spaces(Math.max(i,j)) + '|' + hybridization;
      lower = spaces(Math.max(i,j)-j) + s2_reverse.substring(0,j+1)+lower;
      return;
    } 
    
    // stacking region
    double caseH_1 = e.stacking(i-1, j-1, i, j) + H[i-1][j-1];
    
    // bulge loop (first sequence)
    double caseH_2 = 0.0;  
    int caseH_2_i_prime = 0;
    for (int x=0; x<i-1; x++) {
      double temp = e.bulge(x,j-1,i,j) + H[x][j-1];
      if (temp < caseH_2) {
        caseH_2 = temp;
        caseH_2_i_prime = x;
      }
    }

    // bulge loop (second sequence)
    double caseH_3 = 0.0;
    int caseH_3_j_prime = 0;
    for (int y=0; y<j-1; y++) {
      double temp = e.bulge(i-1,y,i,j)+H[i-1][y];
      if (temp < caseH_3) {
        caseH_3 = temp;
        caseH_3_j_prime = y;
      }
    }
    
    // interior loop
    double caseH_4 = 0.0;
    int caseH_4_i_prime = 0;
    int caseH_4_j_prime = 0;
    for (int x=0; x<i-1; x++) {
      for (int y=0; y<j-1; y++) {
        double temp = e.interior(x,y,i,j) + H[x][y];
        if (temp <caseH_4) {
          caseH_4 = temp;
          caseH_4_i_prime = x;
          caseH_4_j_prime = y;
        }
      }
    }  

 
    if (H[i][j] == caseH_1) { // H(i,j) ends a stacking region
      upper = s1.charAt(i) +upper;
      hybridization = '|'+hybridization;
      lower = s2_reverse.charAt(j)+lower;
      determineHybridization(i-1,j-1);
    }
    
    if (H[i][j] == caseH_2) { // H(i,j) ends a bulge loop in S
      upper = s1.substring(caseH_2_i_prime+1,i+1) + upper;
      hybridization = spaces(i-caseH_2_i_prime-1) + '|' + hybridization;
      lower = gaps(i-caseH_2_i_prime-1) + s2_reverse.charAt(j) + lower;
      determineHybridization(caseH_2_i_prime,j-1);
    } 
    
    if (H[i][j] == caseH_3) { // H(i,j) ends a bulge loop in T
      upper = gaps(j-caseH_3_j_prime-1) + s1.charAt(i)+upper;
      hybridization = spaces(j-caseH_3_j_prime-1) + '|' + hybridization;
      lower = s2_reverse.substring(caseH_3_j_prime+1, j+1) + lower;
      determineHybridization(i-1,caseH_3_j_prime);
    } 
     
    if (H[i][j] == caseH_4) { // H(i,j) ends an interior loop
      int maxLoop = Math.max(i-caseH_4_i_prime-1,j-caseH_4_j_prime-1);
      upper = s1.substring(caseH_4_i_prime+1,i)+gaps(Math.max(maxLoop-(i-caseH_4_i_prime-1),0)) + s1.charAt(i)+upper;
      hybridization = spaces(maxLoop) + '|' + hybridization;
      lower = s2_reverse.substring(caseH_4_j_prime+1,j) + gaps(Math.max(maxLoop-(j-caseH_4_j_prime-1),0))+s2_reverse.charAt(j)+lower;
      determineHybridization(caseH_4_i_prime,caseH_4_j_prime);
    }
  }
  
  // Returns the first sequence (including any gaps) participating in the hybridization
  public String getS1Hybridization() {
    return upper;
  }
  
  // Returns the sequence of vertical bars indicates interacting basepairs from the two sequences
  public String getHybridization() {
    return hybridization;
  }
  
  // Returns the second sequence (including any gaps) participating in the hybridization
  public String getS2Hybridization() {
    return lower;
  }
  
  // private instance methods
  
  // Returns a string of x space characters
  private String spaces(int x) {
    if (x==0) {
      return "";
    }else {
      return String.format("%0"+x+"d", 0).replace("0"," ");
    }
  }
  
  // Returns a string of x gap characters
  private String gaps(int x) {
    if (x==0) {
      return "";
    } else {
      return String.format("%0"+x+"d", 0).replace("0","-");
    }
  }
  
  // Returns minimum of four numbers
  private double min(double a, double b, double c, double d, double e) {
      return Math.min(Math.min(Math.min(a,b), Math.min(c,d)), e);
  }
  
  private String arrayToString() {
    if (H.length == 0) return "";
    if (H[0].length == 0) return "";
    StringBuilder sb = new StringBuilder();
    java.text.DecimalFormat df = new java.text.DecimalFormat("0.0");
    for (int j=s2.length()-1; j>=0; j--)
      sb.append("\t" + s2.charAt(j));
    sb.append("\n");
    for (int i=0; i<s1.length(); i++) {
      sb.append(s1.charAt(i));
      for (int j=0; j<s2.length(); j++) {
        sb.append("\t" + df.format(H[i][j]));
      }
     sb.append("\n");
    }
    
    return sb.toString();
    }
 
  // Main method 
  public static void main(String[] args) {
    if (args.length < 2) {
      System.err.println("\nTo execute this program, two command line arguments are required corresponding to the names of two FASTA files, each containing an RNA sequence. For example,\n");
      System.err.println("\tjava TEST_RNA_Hybridization sequences/RNA_2a.txt sequences/RNA_2b.txt\n");
      return;
    }

    RNA_Hybridization r = new RNA_Hybridization(args[0], args[1]);
    r.computeHybridizationWithMinimumEnergy();
    r.determineHybridization();
    //System.out.println(r.arrayToString());
    System.out.println("\nHybridization with lowest energy:\n");
    System.out.println("\t" + r.getS1Hybridization());
    System.out.println("\t" + r.getHybridization());
    System.out.println("\t" + r.getS2Hybridization());
    System.out.println("\nOptimal hybridization energy: " + r.get_optimalEnergy() + "\n");    
  }  
  }
