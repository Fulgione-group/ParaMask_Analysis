package c.e.data_processing;

import java.io.File;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Scanner;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.*;
import java.io.*;
import java.lang.reflect.Constructor;

public class VcfCombined_to_snpMatrix_alpina {
	
	public VcfCombined_to_snpMatrix_alpina() {}
	
	public void setFileToConvert(String filename, String chromosome, String outFile){
		
		try {
			//
			// Open combined vcf
			File vcf_comb = new File(filename);
			Scanner scannerVcf = new Scanner(vcf_comb);
			
			// Open also the writing file
			Writer writer = new FileWriter(outFile + "_c5q25_chr" + chromosome + ".txt");
			PrintWriter out = new PrintWriter(writer);
			
			int checkChr = 0;
			int checkPos=0;	
			int nSample =0;
			Boolean doublePos = false;
			char[] bases = null;
			int fillBases = 0;
			String[] splitIDs = null;
				
			while ( scannerVcf.hasNextLine() ) {
				String snp1 = scannerVcf.nextLine();
			       	String[] splitSnp = snp1.split("\t");
		       	
				// Mind the header now
				if (snp1.charAt(0) == '#') {
					if (snp1.charAt(1) != '#') {
						splitIDs = snp1.split("\t");
						
						// Count the samples
					       	nSample = splitSnp.length - 9;
					       	System.out.println("Nsamples: " + nSample);
				       		bases = new char[nSample + 1];
			                        for (int bb=0; bb<bases.length; bb++) {
                        			        bases[bb] = 'N';
			                        }
						
				       		// Print out header and samples
				       		if ( (chromosome.charAt(0) == '1' ) ) {
							out.print("chromosome" + "\t" + "position" + "\t" + "ref");
				       			for (int sample=9; sample<splitIDs.length; sample++) {
								out.print( "\t" + splitSnp[sample]);
				       			}
		       					out.print("\n");
						}
					}
				}
				// Header is over
				
				// mind proper stuff now
				
				if (splitSnp[0].charAt(0) == 'c' && splitSnp[0].charAt(1) == 'h' && splitSnp[0].charAt(2) == 'r' && splitSnp[0].charAt(3) == chromosome.charAt(0) ) { 
				
					
						// Indels out - ref length == 1
						if (splitSnp[3].length() == 1) {
							// check also alt length == 1
							String[] check = splitSnp[4].split(",");
							boolean inDel=false;
							for (int c=0; c< check.length; c++) {
								if (check[c].length() > 1) {
									inDel = true;
								}
							}
							
							// Also chlorop and mytoc out
							if (!inDel) {
								doublePos = false;
								// Just check doubles problem
								if ( checkChr == Character.getNumericValue(splitSnp[0].charAt(3))) {						// 3 for our vcf
									if ( (checkPos == Integer.parseInt(splitSnp[1])) && (checkPos != 0) ) {
										//
										// Do something about it
										doublePos = true;
									} else {
										// Out if they are not all Ns
										// If it is the first position, bases[] is all Ns as well
										Boolean allN = true;
										for (int n=1; n<bases.length; n++) {
											if (bases[n] != 'N') {
												allN = false;
											}
										}
										if (!allN) {
											out.print(checkChr + "\t" + checkPos);
                                                                               		for (int bb=0; bb<bases.length; bb++) {
                                                                                        	out.print("\t" + bases[bb]);
                                                                            		}
                                                                                	out.print("\n");
										}
										bases = new char[nSample + 1];
										for (int bb=0; bb<bases.length; bb++) {
                        							        bases[bb] = 'N';
							                        }
									}
									checkPos = Integer.parseInt(splitSnp[1]);
								} else {
									checkChr = Character.getNumericValue(splitSnp[0].charAt(3));						//3 for our vcf
									checkPos = Integer.parseInt(splitSnp[1]);	
								}
							
								// Reference base as bases[0]
								
								bases[0] = splitSnp[3].charAt(0);
								
							       	// For every sample, print a base
							       	int run = 0;
							       	for (int sample=9; sample<splitIDs.length; sample++) {
					       					// Is there any call at all?
					    	 	  			//
					     	  				// Normal vcf: GT:AD:DP:GQ:PGT:PID:PL:PS
										int indGT = 0;
										int indDP = 0;
										int indGQ = 0;
										String[] formatt = splitSnp[8].split(":");
										for (int f=0; f<formatt.length; f++) {
											if (formatt[f].equals("GT")) {
												indGT = f;
											}
											if (formatt[f].equals("DP")) {
												indDP = f;
											}
											if (formatt[f].equals("GQ") || formatt[f].equals("RGQ")) {
												indGQ = f;
											}
										}
										if ( (splitSnp[sample].split(":")[indGT].charAt(0) != '.') && (splitSnp[sample].split(":")[indDP].charAt(0) != '.') && (splitSnp[sample].split(":")[indGQ].charAt(0) != '.') ){
											// Check coverage and quality
							       				int cov=Integer.parseInt(splitSnp[sample].split(":")[indDP]);
						       					int qual=Integer.parseInt(splitSnp[sample].split(":")[indGQ]);
						    		   			
						      			 		// Print out something reasonable
						       					if ( (cov >=5) && (qual >= 25) ) {
								       				// Normal vcf:
								       				if (splitSnp[sample].split(":")[indGT].equals("0/0") || splitSnp[sample].split(":")[indGT].equals("0|0")) {
								       					bases[run + 1] = splitSnp[3].charAt(0);
									       			} else {
									       				// Normal vcf:
									       				if ( (splitSnp[4].length() == 1) && ( splitSnp[sample].split(":")[indGT].equals("1/1") || splitSnp[sample].split(":")[indGT].equals("1|1")) ) {
	                                               		         					bases[run + 1] = splitSnp[4].charAt(0);
									       				} else {
														if ( (splitSnp[4].length() == 1) && (splitSnp[sample].split(":")[indGT].equals("0|1") || splitSnp[sample].split(":")[indGT].equals("1|0") || splitSnp[sample].split(":")[indGT].equals("0/1") || splitSnp[sample].split(":")[indGT].equals("1/0") ) ) {
															// Sorteggia un allele o l'altro per heterozygote sites
															Random randomGenerator = new Random();
															int ran = randomGenerator.nextInt(2);
															bases[run + 1] = splitSnp[3+ran].charAt(0);
														} else {
															bases[run + 1] = 'N';
														}
							       						}
							       					}
											} else {
                                                               					if (!doublePos) {
													bases[run + 1] = 'N';
													// out.print("\t" + "N");
												}
						     	  				}
				       						} else {
											if (!doublePos) {
				       								bases[run + 1] = 'N';
				       								//out.print("\t" + "N");
				       							}
										}
										run = run + 1;
									//}
								}
							       	// out.print("\n");
							  	// Next line!
							}
						}	
					}
       				}	
			// }
			// Put out the last position
			Boolean allN = true;
                        for (int n=1; n<bases.length; n++) {
                                if (bases[n] != 'N') {
                                	allN = false;
                        	}
                        }
			if (!allN) {
                                out.print(checkChr + "\t" + checkPos);
                                for (int bb=0; bb<bases.length; bb++) {
                                	out.print("\t" + bases[bb]);
                                }
                        	out.print("\n");
                        }
			out.close();
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		VcfCombined_to_snpMatrix_alpina vcfCombined_to_snpMatrix_alpina = new VcfCombined_to_snpMatrix_alpina();
		vcfCombined_to_snpMatrix_alpina.setFileToConvert(args[0], args[1], args[2]);
	}
}
