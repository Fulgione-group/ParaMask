/*
 * Copyright (C) 2024 btjeng
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/**
 *
 * @author btjeng
 */

package ParaMask;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;


public class Paramask_Cluster_Seeds {

    public static void main(String[] args) throws FileNotFoundException, IOException {
        System.out.println("Args: " + Arrays.toString(args));
//        Cluster_seeds(new String[]{"--cov", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.stat.txt", "--het", "/home/btjeng/Data/Paramask/testing/Sim0.1HW_rep1EMresults.het", "--cutoff", "322", "--covgw", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.gw.txt"});
//        update_het(new String[]{"--cov", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.stat.txt", "--het", "/home/btjeng/Data/Paramask/testing/Sim0.1HW_rep1EMresults.het", "--cutoff", "322", "--covgw", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.gw.txt"});
//        make_SV_bed(new String[]{"--range", "1,1000000", "--cov", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.stat.txt", "--het", "/home/btjeng/Data/Paramask/testing/Sim0.1HW_rep1EMresults.het", "--cutoff", "322", "--covgw", "/home/btjeng/Data/Paramask/testing/Simulations_0.1_SV_SC_rep1_mat_cpos.vcf.cov.gw.txt"});

        Cluster_seeds(args);
        update_het(args);
        make_SV_bed(args);
    }

    private static void Cluster_seeds(String[] inputList) throws FileNotFoundException, IOException {
        String hetpath = null;
        String covpath = null;
        String covGWpath = null;
        String outpath = null;
        int distcutoff = 0;
        int nSNPpurge = 1;
        for (int i = 0; i < inputList.length; i++) {
            // full input path of the vcf file
            if (inputList[i].equals("--cov") | inputList[i].equals("-c")) {
                covpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--covgw") | inputList[i].equals("-cg")) {
                covGWpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--het") | inputList[i].equals("-h")) {
                hetpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--cutoff") | inputList[i].equals("-cd")) {
                distcutoff = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
            if (inputList[i].equals("--purge") | inputList[i].equals("-p")) {
                nSNPpurge = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
            //full path of the output vcf 
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                outpath = inputList[(i + 1)];
                i++;
            }
        }
        if (outpath == null) {
            outpath = hetpath.substring(0, hetpath.lastIndexOf(".")) + ".clusters.txt";
        }
        InputStream inputStream = new BufferedInputStream(new FileInputStream(covpath));
        BufferedReader readerCov = new BufferedReader(new InputStreamReader(inputStream));
        InputStream inputStream2 = new BufferedInputStream(new FileInputStream(hetpath));
        BufferedReader readerHet = new BufferedReader(new InputStreamReader(inputStream2));
        InputStream inputStream4 = new BufferedInputStream(new FileInputStream(hetpath));
        BufferedReader readerHet2 = new BufferedReader(new InputStreamReader(inputStream4));
        InputStream inputStream3 = new BufferedInputStream(new FileInputStream(covGWpath));
        BufferedReader readerCovGW = new BufferedReader(new InputStreamReader(inputStream3));
        BufferedWriter writer = null;
        writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))));
        ArrayList<int[]> clusters = new ArrayList<>();
        int linecounter = 0;
        String[] linearray = null;
        Float[] covGW = null;
        String lineCovGW = null;
        while ((lineCovGW = readerCovGW.readLine()) != null) {
            if (linecounter == 1) {
                covGW = Arrays.stream(P_TAB.split(lineCovGW)).map(Float::valueOf).toArray(Float[]::new);
            }
            linecounter++;
        }
//        System.out.println(Arrays.toString(covGW));
        String lineHet2 = null;
        int nrow = 0;
        linecounter = 0;
        while ((lineHet2 = readerHet2.readLine()) != null) {
            linecounter++;
        }
        nrow = linecounter;
        boolean extend = false;
        boolean extendl = false;
        boolean bridge = false;
        int ccounter = 1;
        // Assuming set1 is a list of objects with properties Position, EM_class, and Mean.coverage
        // You need to replace the type of set1 with the actual type in your code.
        // Assuming nrow is a method that returns the number of rows in set1
        linecounter = 0;
        String lineHet = null;
        String chrom = null;
        String lineCov = null;
        String prevsubsample = null;
        String subsample = null;
        String[] lineCovarray = null;
        String[] lineHetarray = null;
        String[] prevlineCovarray = null;
        String[] prevlineHetarray = null;
        int pos = 0;
        int prevpos = 0;
        int posafter = 0;
        float tmpCov = 0;
        int nonnaSubSet = 0;
        String[] subsetCovArray = null;
        String[] subsampleArray = null;
        String[] prevsubsampleArray = null;
        float subsetCovGW = 0;
        float subsetCov = 0;
        String[] tmpBridge = null;
        //still need to assign the lists after each cluster has been broken or at the beginning until the first seed
        ArrayList<String[]> listCovarray = new ArrayList<>();
        ArrayList<String[]> listHetarray = new ArrayList<>();
        HashMap<String, Integer> samplesIndex = new HashMap<String, Integer>();
        ArrayList<String[]> tmpCluster = null;
        writer.write("chromosome" + "\t" + "position" + "\t" + "cluster" + "\t" + "emClass" + "\t" + "clusterCause" + "\t" + "CovHetOfLastSeed" + "\t" + "CovGWHetOfLastSeed" + "\t" + "hetGenOfLastSeed");
        writer.newLine();
        while ((lineHet = readerHet.readLine()) != null) {
//            System.out.println(subsetCov); 
            lineCov = readerCov.readLine();
            if (linecounter > 0) {
                lineCovarray = P_TAB.split(lineCov);
                lineHetarray = P_TAB.split(lineHet);
//                System.out.println(Arrays.toString(lineCovarray));
                if (!lineCovarray[1].equals(lineHetarray[1])) {
                    throw new IllegalArgumentException("positions do not match. Check .Het and .Cov files");
                }
                pos = Integer.valueOf(lineCovarray[1]);
                chrom = lineCovarray[0];
                if (linecounter > 1) {
                    prevpos = Integer.valueOf(prevlineCovarray[1]);
                }
                int emClass = Integer.valueOf(lineHetarray[17]);
                int bridgeDist = 0;
                if (!extend) {
                    if (emClass == 2) {
//                        ccounter++;
                        prevsubsample = lineCovarray[(lineCovarray.length - 1)].replaceAll("\\,", ":");
                        prevsubsampleArray = P_COMMA.split(lineCovarray[(lineCovarray.length - 1)]);
                        nonnaSubSet = 0;
                        for (int i = 0; i < prevsubsampleArray.length; i++) {
//                                System.out.println(samplesIndex.get(prevsubsampleArray[i]));
                            tmpCov = Float.parseFloat(lineCovarray[samplesIndex.get(prevsubsampleArray[i])]);
                            if (tmpCov >= 0) {
                                nonnaSubSet++;
                                subsetCovGW += covGW[(samplesIndex.get(prevsubsampleArray[i]) - 2)];
                                subsetCov += Float.parseFloat(lineCovarray[samplesIndex.get(prevsubsampleArray[i])]);
                            }
                        }
                        subsetCovGW = subsetCovGW / nonnaSubSet;
                        subsetCov = subsetCov / nonnaSubSet;
                        extend = true;
                        extendl = true;
                        bridge = false;
                        int j = listCovarray.size();
                        tmpCluster = new ArrayList<>();
                        tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                        subsetCovGW = 0;
                        subsetCov = 0;
                        while (extendl) {
                            j--;
                            if (j < 0) {
                                extendl = false;
                                break;
                            }
                            posafter = pos;
                            pos = Integer.valueOf(listCovarray.get(j)[1]);
                            chrom = listCovarray.get(j)[0];
                            emClass = Integer.valueOf(listHetarray.get(j)[17]);
                            subsetCovArray = listCovarray.get(j);
                            subsampleArray = P_COMMA.split(lineCovarray[(lineCovarray.length - 1)]);
                            subsample = lineCovarray[(lineCovarray.length - 1)].replaceAll("\\,", ":");
                            //need to implement a NA for coverage if the sample is missing!!!!!!!! TBE
                            tmpCov = 0;
                            subsetCov = 0;
                            subsetCovGW = 0;
                            nonnaSubSet = 0;
                            for (int i = 0; i < prevsubsampleArray.length; i++) {
//                                System.out.println(samplesIndex.get(prevsubsampleArray[i]));
                                tmpCov = Float.parseFloat(subsetCovArray[samplesIndex.get(prevsubsampleArray[i])]);
                                if (tmpCov >= 0) {
                                    nonnaSubSet++;
                                    subsetCovGW += covGW[(samplesIndex.get(prevsubsampleArray[i]) - 2)];
                                    subsetCov += Float.parseFloat(subsetCovArray[samplesIndex.get(prevsubsampleArray[i])]);
                                }
                            }
                            //subsetCov and subsetCovGW become NA when all genotypes are missing and subsetCov >= (1.5 * subsetCovGW) becomes false
                            subsetCovGW = subsetCovGW / nonnaSubSet;
                            subsetCov = subsetCov / nonnaSubSet;
//                            System.out.println(subsetCov); 
//                            if (nonnaSubSet == 0) {
//                                System.out.println(subsetCov >= (1.5 * subsetCovGW));
//                            }
                            if (!bridge) {
                                if (emClass == 2) {
                                    prevsubsampleArray = subsampleArray;
                                    prevsubsample = subsample;
                                    tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                } else if (subsetCov >= (1.5 * subsetCovGW)) {
                                    tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "hetcov", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                } else if (emClass == 1) {
                                    bridge = true;
                                    bridgeDist = posafter - pos;
                                    tmpBridge = new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "bridge", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample};
                                } else {
                                    extendl = false;
                                    break;
                                }
                            } else {
                                int tmpDist = posafter - pos;
                                if ((float) (tmpDist + bridgeDist) / 2 < distcutoff) {
                                    if (emClass == 2) {
                                        prevsubsampleArray = subsampleArray;
                                        prevsubsample = subsample;
                                        tmpCluster.add(tmpBridge);
                                        tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                        bridge = false;
                                    } else if (subsetCov >= (1.5 * subsetCovGW)) {
                                        tmpCluster.add(tmpBridge);
                                        tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "hetcov", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                        bridge = false;
//                                    } else if (emClass == 1) {
//                                        tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), "bridge", prevsubsample});
//                                        bridge = false;
                                    } else {
                                        extendl = false;
                                        bridge = false;
                                        break;
                                    }
                                } else {
                                    extendl = false;
                                    bridge = false;
                                    break;
                                }
                            }
                        }

                        int maxRow = tmpCluster.size();
                        if (!tmpCluster.isEmpty()) {
                            // Reverse order
                            Collections.reverse(tmpCluster);
                        }
                        bridge = false;
                    }
                    listCovarray.add(lineCovarray);
                    listHetarray.add(lineHetarray);
                } else {
                    subsampleArray = P_COMMA.split(lineCovarray[(lineCovarray.length - 1)]);
                    subsample = lineCovarray[(lineCovarray.length - 1)].replaceAll("\\,", ":");;
                    tmpCov = 0;
                    subsetCov = 0;
                    subsetCovGW = 0;
                    nonnaSubSet = 0;
//                    System.out.println(subsetCov);
                    for (int i = 0; i < prevsubsampleArray.length; i++) {
//                                System.out.println(samplesIndex.get(prevsubsampleArray[i]));
                        tmpCov = Float.parseFloat(lineCovarray[samplesIndex.get(prevsubsampleArray[i])]);
                        if (tmpCov >= 0) {
                            nonnaSubSet++;
                            subsetCovGW += covGW[(samplesIndex.get(prevsubsampleArray[i]) - 2)];
                            subsetCov += Float.parseFloat(lineCovarray[samplesIndex.get(prevsubsampleArray[i])]);
//                            System.out.println(covGW[(samplesIndex.get(prevsubsampleArray[i]) - 2)]);
//                            System.out.println(Arrays.toString(lineCovarray));
                        }
                    }
//                    System.out.println(subsetCov);
//                    System.out.println(subsetCovGW);
//                    System.out.println(nonnaSubSet);
                    subsetCovGW = subsetCovGW / nonnaSubSet;
                    subsetCov = subsetCov / nonnaSubSet;
//                    System.out.println(subsetCov); 
//                    System.out.println(subsetCovGW); 
                    if (!bridge) {
                        if (emClass == 2) {
                            prevsubsampleArray = subsampleArray;
                            prevsubsample = subsample;
                            tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                            if (linecounter == (nrow - 1)) {
                                if (tmpCluster.size() > nSNPpurge) {
                                    for (int t = 0; t < tmpCluster.size(); t++) {
                                        writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                        writer.newLine();
                                    }
                                    ccounter++;
                                }
                                tmpCluster.clear();
                            }
                        } else if (subsetCov >= (1.5 * subsetCovGW)) {
                            tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "hetcov", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                            if (linecounter == (nrow - 1)) {
                                if (tmpCluster.size() > nSNPpurge) {
                                    for (int t = 0; t < tmpCluster.size(); t++) {
                                        writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                        writer.newLine();
                                    }
                                    ccounter++;
                                }
                                tmpCluster.clear();
                            }
                            bridge = false;
                        } else if (emClass == 1) {
                            bridge = true;
                            bridgeDist = pos - prevpos;
                            tmpBridge = new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "bridge", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample};
                            if (linecounter == (nrow - 1)) {
                                if (tmpCluster.size() > nSNPpurge) {
                                    for (int t = 0; t < tmpCluster.size(); t++) {
                                        writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                        writer.newLine();
                                    }
                                    ccounter++;
                                }
                                tmpCluster.clear();
                            }
                        } else {
                            if (tmpCluster.size() > nSNPpurge) {
                                for (int t = 0; t < tmpCluster.size(); t++) {
                                    writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                    writer.newLine();
                                }
                                ccounter++;
                            }
                            extend = false;
                            listCovarray = new ArrayList<>();
                            listHetarray = new ArrayList<>();
                            // break();
                        }
                    } else {
                        int tmpDist = pos - prevpos;
                        if ((tmpDist + bridgeDist) / 2 < distcutoff) {
                            if (emClass == 2) {
                                prevsubsampleArray = subsampleArray;
                                prevsubsample = subsample;
                                tmpCluster.add(tmpBridge);
                                tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                bridge = false;
                                if (linecounter == (nrow - 1)) {
                                    if (tmpCluster.size() > nSNPpurge) {
                                        for (int t = 0; t < tmpCluster.size(); t++) {
                                            writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                            writer.newLine();
                                        }
                                        ccounter++;
                                    }
                                    tmpCluster.clear();
                                }
                            } else if (subsetCov >= (1.5 * subsetCovGW)) {
                                tmpCluster.add(tmpBridge);
                                tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "hetcov", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                bridge = false;
                                if (linecounter == (nrow - 1)) {
                                    if (tmpCluster.size() > nSNPpurge) {
                                        for (int t = 0; t < tmpCluster.size(); t++) {
                                            writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                            writer.newLine();
                                        }
                                        ccounter++;
                                    }
                                    tmpCluster.clear();
                                }
                            } else {
                                if (tmpCluster.size() > nSNPpurge) {
                                    for (int t = 0; t < tmpCluster.size(); t++) {
                                        writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                        writer.newLine();
                                    }
                                    ccounter++;
                                }
                                extend = false;
                                bridge = false;
                                listCovarray = new ArrayList<>();
                                listCovarray.add(lineCovarray);
                                listHetarray = new ArrayList<>();
                                listHetarray.add(lineHetarray);
                                // break();
                            }
                        } else {
                            if (tmpCluster.size() > nSNPpurge) {
                                for (int t = 0; t < tmpCluster.size(); t++) {
                                    writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                    writer.newLine();
                                }
                                ccounter++;
                            }
                            if (emClass == 2) {
//                                ccounter++;
                                tmpCluster = new ArrayList<>();
                                prevsubsampleArray = subsampleArray;
                                prevsubsample = subsample;
                                tmpCluster.add(new String[]{chrom, Integer.toString(pos), Integer.toString(ccounter), String.valueOf(emClass), "seed", String.valueOf(subsetCov), String.valueOf(subsetCovGW), prevsubsample});
                                bridge = false;
                                extend = true;
                                listCovarray = new ArrayList<>();
                                listHetarray = new ArrayList<>();
                                if (linecounter == (nrow - 1)) {
                                    if (tmpCluster.size() > nSNPpurge) {
                                        writer.write(Arrays.toString(tmpCluster.get(0)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                                        writer.newLine();
                                        ccounter++;
                                    }
                                    tmpCluster.clear();
                                }
                            } else {
                                extend = false;
                                bridge = false;
                                listCovarray = new ArrayList<>();
                                listCovarray.add(lineCovarray);
                                listHetarray = new ArrayList<>();
                                listHetarray.add(lineHetarray);
                            }
                        }
                    }
                }
                prevlineCovarray = lineCovarray;
                prevlineHetarray = lineHetarray;
            } else {
                lineCovarray = P_TAB.split(lineCov);
                for (int i = 2; i < (lineCovarray.length - 1); i++) {
                    samplesIndex.put(lineCovarray[i], i);
//                    System.out.println(lineCovarray[i] + "\t" + i);
                }
            }
            linecounter++;
        }
        if (extend) {
            // Add the last cluster if needed
            if (tmpCluster.size() > nSNPpurge) {
                for (int t = 0; t < tmpCluster.size(); t++) {
                    writer.write(Arrays.toString(tmpCluster.get(t)).replaceAll("\\[|\\]|\\p{Zs}+", "").replaceAll("\\,", "\t"));
                    writer.newLine();
                }
                ccounter++;
            }
        }
        writer.flush();
        writer.close();
    }

    private static void update_het(String[] inputList) throws FileNotFoundException, IOException {
        String hetpath = null;
        String hetpathnew = null;
        String inpath = null;
        for (int i = 0; i < inputList.length; i++) {
            if (inputList[i].equals("--het") | inputList[i].equals("-h")) {
                hetpath = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--out") | inputList[i].equals("-o")) {
                inpath = inputList[(i + 1)];
                i++;
            }
        }
        if (inpath == null) {
            inpath = hetpath.substring(0, hetpath.lastIndexOf(".")) + ".clusters.txt";
        }
        InputStream inputStream4 = new BufferedInputStream(new FileInputStream(hetpath));
        BufferedReader readerHet = new BufferedReader(new InputStreamReader(inputStream4));
        InputStream inputStream = new BufferedInputStream(new FileInputStream(inpath));
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        BufferedWriter writer = null;
        hetpathnew = hetpath.substring(0, hetpath.lastIndexOf(".")) + ".finalClass.het";
        writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(hetpathnew))));
        String line = null;
        String[] linearray = null;
        ArrayList<Integer> parPos = new ArrayList<>();
        ArrayList<Integer> parCluster = new ArrayList<>();
        int linecounter = 0;
        while ((line = reader.readLine()) != null) {
            if (linecounter > 0) {
                linearray = P_TAB.split(line);
                parPos.add(Integer.valueOf(linearray[1]));
                parCluster.add(Integer.valueOf(linearray[2]));
//                System.out.println(parPos.get((linecounter-1)));
            }
            linecounter++;
        }
        int parPoscounter = 0;
        linecounter = 0;
        while ((line = readerHet.readLine()) != null) {
            if (linecounter > 0) {
                linearray = P_TAB.split(line);
                boolean ispos = parPoscounter != parPos.size() && Integer.valueOf(linearray[1]).equals(parPos.get(parPoscounter));
//                System.out.println(linearray[1] + "\t" + parPos.get(parPoscounter) + "\t" + ispos);
                if (ispos) {
                    writer.write(line + "\t" + "1" + "\t" + parCluster.get(parPoscounter));
                    writer.newLine();
                    parPoscounter++;
                } else {
                    writer.write(line + "\t" + "0" + "\t" + "0");
                    writer.newLine();
                }
            } else {
                writer.write(line + "\t" + "finalClass" + "\t" + "cluster");
                writer.newLine();
            }
            linecounter++;
        }
        writer.flush();
        writer.close();
    }

    private static void make_SV_bed(String[] inputList) throws FileNotFoundException, IOException {
        String hetpath = null;
        String hetpathnew = null;
        String inpath = null;
        String range = null;
        int distcutoff = 0;
        for (int i = 0; i < inputList.length; i++) {
            // full input path of the vcf file
            if (inputList[i].equals("--het") | inputList[i].equals("-h")) {
                hetpath = inputList[(i + 1)];
                i++;
            }
            //full path of the output vcf 
            if (inputList[i].equals("--range") | inputList[i].equals("-r")) {
                range = inputList[(i + 1)];
                i++;
            }
            if (inputList[i].equals("--cutoff") | inputList[i].equals("-c")) {
                distcutoff = Integer.parseInt(inputList[(i + 1)]);
                i++;
            }
        }
        hetpathnew = hetpath.substring(0, hetpath.lastIndexOf(".")) + ".finalClass.het";
        Integer chrstart = 0;
        Integer chrend = 0;
        String[] rangeArray = P_COMMA.split(range);
        chrstart = Integer.valueOf(rangeArray[0]);
        chrend = Integer.valueOf(rangeArray[1]);
        String outpath = hetpathnew.substring(0, hetpathnew.lastIndexOf(".")) + ".bed";
        InputStream inputStream = new BufferedInputStream(new FileInputStream(hetpathnew));
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        BufferedWriter writer = null;
        writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new java.io.File(outpath))));
        writer.write("Chromosome" + "\t" + "Start" + "\t" + "End" + "\t" + "type:0-single-copy;1-multi-copy" + "\t" + "nSNPs" + "\t" + "Cluster");
        writer.newLine();
        int start = 0;
        int end = 0;
        int nSNPs = 0;
        ArrayList<String[]> SVTabInfer = new ArrayList<>();
        String line = null;
        String chrom = null;
        String[] linearray = null;
        // Assuming set1 is an object of a class with the necessary properties
        int linecounter = 0;
        int distcutoffhalv = (int) Math.floor(distcutoff / 2);
        int distextend = 10000;
        String state = null;
        String cluster = null;
        while ((line = reader.readLine()) != null) {
            linearray = P_TAB.split(line);
            if (linecounter > 0) {
                if (linecounter == 1) {
                    state = linearray[19];
                    chrom = linearray[0];
                    start = Integer.valueOf(linearray[1]);
                    start = chrstart;
                } else if (!state.equals(linearray[19])) {
                    int distSNP = Integer.valueOf(linearray[1]) - end;
                    distextend = Math.min(distcutoff, distSNP);
                    int distextendhalv = (int) Math.floor(Float.valueOf(distextend) / 2);
                    if (distextend >= distSNP) {
                        if ((distextend % 2) != 0) {
                            end = end + distextendhalv;
                        } else {
                            end = end + distextendhalv - 1;
                        }
                        writer.write(chrom + "\t" + start + "\t" + end + "\t" + state + "\t" + nSNPs + "\t" + cluster);
                        writer.newLine();
                        chrom = linearray[0];
                        start = Integer.valueOf(linearray[1]);
                        start = start - distextendhalv;
                    } else {
                        if (state.equals("0")) {
                            end = Integer.valueOf(linearray[1]) - distextendhalv - 1;
                            writer.write(chrom + "\t" + start + "\t" + end + "\t" + state + "\t" + nSNPs + "\t" + cluster);
                            writer.newLine();
                            start = end + 1;
                        } else {
                            if ((distextend % 2) != 0) {
                                end = end + distextendhalv;
                            } else {
                                end = end + distextendhalv - 1;
                            }
                            writer.write(chrom + "\t" + start + "\t" + end + "\t" + state + "\t" + nSNPs + "\t" + cluster);
                            writer.newLine();
                            start = end + 1;
                        }

                    }
                    nSNPs = 0;
                }
                cluster = linearray[20];
                nSNPs++;
                state = linearray[19];
                end = Integer.valueOf(linearray[1]);
            }
            linecounter++;
        }
        end = chrend;
        writer.write(chrom + "\t" + start + "\t" + end + "\t" + state + "\t" + nSNPs + "\t" + cluster);
        // Output
        for (String[] row : SVTabInfer) {
            System.out.println(String.join(", ", row));
        }
        writer.flush();
        writer.close();

    }
    private static final Pattern P_TAB = Pattern.compile("\t");
    private static final Pattern P_COMMA = Pattern.compile(",");
}
