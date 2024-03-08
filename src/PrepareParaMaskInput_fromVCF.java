/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Evolutionary_genetics;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author btjeng
 */
public class PrepareParaMaskInput_fromVCF {

    public static void main(String[] args) {
        System.out.println("Args: " + Arrays.toString(args));
//        Afreq_vs_hets(new String[]{"--vcf", "/home/btjeng/Data/testing/test.vcf", "-VF"});
// to build java byte code
        Afreq_vs_hets(args);
    }

    private static void Afreq_vs_hets(String[] inputList) {
        //start with input parameters
        String vcfpath1 = null;
        String outpath = null;
        String outpath2 = null;
        String outpath3 = null;
        String popfile = null;
        float missingness = 0;
        BufferedReader readerPop = null;
        boolean VF = true;
        for (int i = 0; i < inputList.length; i++) {
            //full input path for the vcf file
            if (inputList[i].equals("--vcf") | inputList[i].equals("-v")) {
                vcfpath1 = inputList[(i + 1)];
                i++;
            }
            // full path for the output
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)] + ".het.stat.txt";
                outpath2 = inputList[(i + 1)] + ".cov.stat.txt";
                outpath3 = inputList[(i + 1)] + ".cov.gw.txt";
                i++;
            }
            // full path to a popfile. If specified every individual in this file will be analyzed. Sample IDs are writen on seperate lines
            if (inputList[i].equals("--popfile") | inputList[i].equals("-p")) {
                popfile = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--missingness") | inputList[i].equals("-m")) {
                missingness = Float.parseFloat(inputList[(i + 1)]);
                i++;
            }
            // if specified the format is assumed to vary across positions
            if (inputList[i].equals("--noVaryingFormat") | inputList[i].equals("-nVF")) {
                VF = false;
            }
        }
        // if no output path is provided it will be constructed from the input file
        if (outpath == null) {
            if (vcfpath1.endsWith(".gz")) {
                outpath = vcfpath1.substring(0, vcfpath1.lastIndexOf(".")) + ".het.stat.txt";
                outpath2 = vcfpath1.substring(0, vcfpath1.lastIndexOf(".")) + ".het.cov.stat.txt";
                outpath3 = vcfpath1.substring(0, vcfpath1.lastIndexOf(".")) + ".cov.gw.txt";
            } else {
                outpath = vcfpath1 + ".het.stat.txt";
                outpath2 = vcfpath1 + ".cov.stat.txt";
                outpath3 = vcfpath1 + ".cov.gw.txt";
            }
        }
        if (vcfpath1 == null) {
            throw new IllegalArgumentException("VCF not specified");
        }
        try {
            // create in- and output streams
            BufferedWriter writer = null;
            BufferedWriter writer2 = null;
            BufferedWriter writer3 = null;
            writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))));
            writer2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath2))));
            writer3 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath3))));
            InputStream inputStream = new BufferedInputStream(new FileInputStream(vcfpath1));
            if (vcfpath1.endsWith(".gz")) {
                inputStream = new GZIPInputStream(inputStream);
            }
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            if (popfile != null) {
                readerPop = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(popfile))));
            }
            List<String> samples = new ArrayList<>();
            List<Integer> samplepos = new ArrayList<>();
            ArrayList<Integer[]> sampleHist = new ArrayList<Integer[]>();
            float Het = 0;
            float Hom1 = 0;
            float Hom2 = 0;
            float Afreq = 0;
            int Nonmissing = 0;
            int refADhet = 0;
            int altADhet = 0;
            int totalhetdepth = 0;
            float allelicRatio = 0;
            float allelicDeviation = 0;
            float meancov = 0;
            float meancovhet = 0;
            float meancovhom = 0;
            String[] linearray = null;
            int[] covarray = null;
            double[] covarraygw = null;
            double[] nsites = null;
            String line;
            List<List<Integer>> popposBin = new ArrayList<>();
            List<String> hetlist = new ArrayList<>();
            //int should be sufficient
//            System.out.println(noOfcomp);
            //haploid sample size
            String chromosome = "start";
            boolean altfound = false;
            int altpos = 4;
            int linecounter = 0;
            List<String> poplist = new ArrayList<>();
            //read the popfile
            if (popfile != null) {
                while ((line = readerPop.readLine()) != null) {
                    poplist.add(line);
                }
            }
            boolean formatfound = false;
            boolean formatcodefound = false;
            boolean formatdpfound = false;
            boolean formatadfound = false;
            int formatcode = 0;
            int adindex = 0;
            int formatpos = 0;
            String[] formatarray = null;
            reader.mark(100000000);
            //search for the format coding. If --VaryingFormat is not specified it will be assumed the format is the same at every position in the VCF
            while (!formatfound) {
                line = reader.readLine();
                linecounter++;
                if (line.startsWith("#CHROM")) {
                    linearray = P_TAB.split(line);
                    for (int k = 0; k < linearray.length; k++) {
                        if (linearray[k].equals("FORMAT")) {
                            formatpos = k;
                            formatcodefound = true;
                        }
                    }
                } else if (formatcodefound) {
                    linearray = P_TAB.split(line);
                    formatarray = linearray[formatpos].split(":");
                    for (int k = 0; k < formatarray.length; k++) {
                        if (formatarray[k].equals("DP")) {
                            formatcode = k;
                            formatdpfound = true;
                        } else if (formatarray[k].equals("AD")) {
                            adindex = k;
                            formatadfound = true;
                        }
                        if (formatdpfound & formatadfound) {
                            formatfound = true;
                            break;
                        }
                    }

                }
            }
            reader.reset();
            int offset = (formatpos + 1);
            linecounter = 0;
            //write the header of the output file
            writer.write("Chromosome" + "\t" + "Position" + "\t" + "Non-missing" + "\t" + "Minor-allele-freq" + "\t" + "Heterozygous-geno-freq" + "\t" + "Homozygous1-geno-freq" + "\t" + "Homozygous2-geno-freq" + "\t" + "Mean-coverage" + "\t" + "Mean-coverage-hom" + "\t" + "Mean-coverage-het" + "\t" + "Het-reference-allele-depth" + "\t" + "Het-alt-allele-depth" + "\t" + "Het-allele-ratio" + "\t" + "Het-allele-deviation");
            writer.newLine();
            //go through the VCF
            while ((line = reader.readLine()) != null) {
                linecounter++;
//                System.out.println(linecounter);
                List<String> pop1tmp = new ArrayList<>();
                //analyze only non header lines.
                if (line.charAt(0) != '#') {
                    linearray = P_TAB.split(line);
                    //make sure to exclude invariant sites or multi-allelic sites
                    if (!linearray[altpos].equals(".") & !linearray[altpos].contains(",")) {
                        if (VF == true) {
                            formatarray = linearray[formatpos].split(":");
                            for (int k = 0; k < formatarray.length; k++) {
                                if (formatarray[k].equals("DP")) {
                                    formatcode = k;
                                }
                                if (formatarray[k].equals("AD")) {
                                    adindex = k;
                                }
                            }
                        }
                        int newstart = 0;
                        if (popfile != null) {
                            //go through every genotype at a certain position of samples specified in the popfile
                            for (int i = 0; i < samplepos.size(); i++) {
                                if (!linearray[samplepos.get(i)].startsWith(".")) {
                                    Nonmissing++;
                                    // retrieve the alleles of the two chromosomes
                                    int g1 = Integer.parseInt(linearray[samplepos.get(i)].split("/|\\||:")[0]);
                                    int g2 = Integer.parseInt(linearray[samplepos.get(i)].split("/|\\||:")[1]);
                                    Afreq += g1 + g2;
                                    meancov += Integer.parseInt(linearray[samplepos.get(i)].split(":")[formatcode]);
                                    covarray[i] = Integer.parseInt(linearray[samplepos.get(i)].split(":")[formatcode]);
                                    covarraygw[i] += Integer.parseInt(linearray[samplepos.get(i)].split(":")[formatcode]);
                                    nsites[i]++;
//                                    System.out.println(covarray[i]);
                                    if (g1 != g2) {
                                        //if the two chromosomes are different count it as heterozygous
                                        Het++;
                                        hetlist.add(samples.get(i));
                                        // also calculate mean coverage (of all heterozygous samples) at that position 
                                        meancovhet += Integer.parseInt(linearray[samplepos.get(i)].split(":")[formatcode]);
                                        if (g1 == 0) {
                                            refADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[0]);
                                            altADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[1]);
                                        } else {
                                            refADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[1]);
                                            altADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[0]);
                                        }
                                    } else {
                                        // calculate mean coverage of all homozygous samples at that position 
                                        meancovhom += Integer.parseInt(linearray[samplepos.get(i)].split(":")[formatcode]);
                                        if (g1 == 0) {
                                            Hom1++;
                                        } else {
                                            Hom2++;
                                        }
                                    }
                                }
                            }
                        } else {
                            // go through all samples if popfile is not specified
                            for (int i = offset; i < linearray.length; i++) {
                                if (!linearray[i].startsWith(".")) {
                                    Nonmissing++;
                                    int g1 = Integer.parseInt(linearray[i].split("/|\\||:")[0]);
                                    int g2 = Integer.parseInt(linearray[i].split("/|\\||:")[1]);
                                    Afreq += g1 + g2;
//                                System.out.println(g1 + "\t" + g2 + "\t"+ linearray[i] + "\t" +linearray[8]);
                                    meancov += Integer.parseInt(linearray[i].split(":")[formatcode]);
                                    covarray[(i - offset)] = Integer.parseInt(linearray[i].split(":")[formatcode]);
                                    covarraygw[(i - offset)] += Integer.parseInt(linearray[i].split(":")[formatcode]);
                                    nsites[(i - offset)]++;
//                                    System.out.println(covarray[(i-offset)]);
                                    if (g1 != g2) {
                                        Het++;
                                        hetlist.add(samples.get(i - offset));
                                        meancovhet += Integer.parseInt(linearray[i].split(":")[formatcode]);
                                        if (g1 == 0) {
                                            refADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[0]);
                                            altADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[1]);
                                        } else {
                                            refADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[1]);
                                            altADhet += Integer.parseInt(linearray[i].split(":")[adindex].split(",")[0]);
                                        }
                                    } else {
                                        meancovhom += Integer.parseInt(linearray[i].split(":")[formatcode]);
                                        if (g1 == 0) {
                                            Hom1++;
                                        } else {
                                            Hom2++;
                                        }
                                    }
                                }
                            }
                        }
//                        System.out.println(meancov);
                        //calculate allele frequency at a given position
                        Afreq = Afreq / (2 * Nonmissing);
                        if (Afreq > 0.5) {
                            Afreq = 1 - Afreq;
                        }
                        if ((Nonmissing - Het) > 0) {
                            meancovhom = meancovhom / (Nonmissing - Het);
                        }
                        meancovhet = meancovhet / Het;
                        // just in case the total depth is different from the total depth in AD column like sometimes GATKS
                        totalhetdepth = refADhet + altADhet;
                        allelicRatio = (float) refADhet / totalhetdepth;
                        allelicDeviation = (float) (((float) totalhetdepth / 2 - altADhet) / Math.sqrt(0.25 * (double) totalhetdepth));
                        Het = Het / Nonmissing;
                        Hom1 = Hom1 / Nonmissing;
                        Hom2 = Hom2 / Nonmissing;
                        meancov = meancov / Nonmissing;
                        if (Afreq > 0 & Nonmissing >= samples.size() * (1 - missingness)) {
                            writer.write(linearray[0] + "\t" + linearray[1] + "\t" + Nonmissing + "\t" + Afreq + "\t" + Het + "\t" + Hom1 + "\t" + Hom2 + "\t" + meancov + "\t" + meancovhom + "\t" + meancovhet + "\t" + refADhet + "\t" + altADhet + "\t" + allelicRatio + "\t" + allelicDeviation);
                            writer.newLine();
                            writer2.write(linearray[0] + "\t" + linearray[1] + "\t" + Arrays.toString(covarray).replaceAll("\\[|\\]", "").replaceAll(",", "\t") + "\t" + hetlist.toString().replaceAll("\\[|\\]", "").replaceAll("\\ ", ""));
                            writer2.newLine();
//                            System.out.println(Arrays.toString(covarray));
                        }
                        //reset all stats for next het SNP
                        Hom1 = 0;
                        Hom2 = 0;
                        Het = 0;
                        Afreq = 0;
                        Nonmissing = 0;
                        meancovhet = 0;
                        meancovhom = 0;
                        meancov = 0;
                        refADhet = 0;
                        altADhet = 0;
                        allelicRatio = 0;
                        allelicDeviation = 0;
                        totalhetdepth = 0;
                        covarray = new int[samples.size()];
                        Arrays.fill(covarray, -1);
                        hetlist = new ArrayList<>();
//                        System.out.println(covarray.length);
                    }
                }
                //compare sample names to header line with sample IDs and assign column positions for each sample (in the popfile)
                if (line.startsWith("#CHROM")) {
                    linearray = P_TAB.split(line);
                    for (int k = 0; k < linearray.length; k++) {
                        linearray = P_TAB.split(line);
                        if (linearray[k].equals("ALT")) {
                            altpos = k;
                        }
                        if (linearray[k].equals("FORMAT")) {
                            formatpos = k;
                        }
                    }
                    for (int k = 0; k < linearray.length; k++) {
                        if (k > formatpos) {
                            if (popfile != null) {
                                for (int i = 0; i < poplist.size(); i++) {
                                    if (poplist.get(i).equals(linearray[k])) {
                                        samples.add(linearray[k]);
                                        samplepos.add(k);
//                                        System.out.println(linearray[k]);
                                    }
                                }
                            } else {
                                samples.add(linearray[k]);
                            }
                        }
                    }
                    //write header of cov file
                    writer2.write("Chromosome" + "\t" + "Position" + "\t" + String.join("\t", samples) + "\t" + "Het-individuals");
                    writer2.newLine();
                    covarray = new int[samples.size()];
                    Arrays.fill(covarray, -1);
                    writer3.write(String.join("\t", samples));
                    writer3.newLine();
                    covarray = new int[samples.size()];
                    Arrays.fill(covarray, -1);
                    covarraygw = new double[samples.size()];
                    nsites = new double[samples.size()];
                }
            }
            for (int i = 0; i < covarraygw.length; i++) {
                covarraygw[i] = covarraygw[i] / nsites[i];
            }
            writer3.write(Arrays.toString(covarraygw).replaceAll("\\[|\\]", "").replaceAll(",", "\t"));
            writer3.newLine();
            writer.flush();
            writer.close();
            writer2.flush();
            writer2.close();
            writer3.flush();
            writer3.close();
        } catch (IOException ex) {
            Logger.getLogger(Main.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static final Pattern P_TAB = Pattern.compile("\t");

}
