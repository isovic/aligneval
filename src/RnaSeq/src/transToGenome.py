import utility_sam
import copy
import sys
import transExtract
import re;

def checkBitMask(name, bitMask, tid_regions):
    bmLen = sum([1 for c in bitMask if c == '1' or c == '0'])
    if bmLen == len(bitMask) and bmLen == len(tid_regions[name]):
        return True
    return False

def extractRegions2(tid, tid_regions):
    fields = tid.split("_")
    if len(fields) != 1:
        name = "_".join(fields[0:len(fields)-1])
        bm = fields[-1]
        if checkBitMask(name, bm, tid_regions):
            return name, [region for i, region in enumerate(tid_regions[name]) if bm[i] == '1']
    return tid, tid_regions[tid]

def extractRegions(tid, tid_regions):
    if "_" not in tid:
        return tid, tid_regions[tid]
    name,bm = tid.split("_")
    return name, [region for i, region in enumerate(tid_regions[name]) if bm[i] == '1']

def solve(sam_lines, transToSeq, tid_regions, tid_exons):
    newSam = []
    for line in sam_lines:
        if line.rname == "*":
            newSam.append(line)
            continue

        name, regions = extractRegions(line.rname, tid_regions)
        exons = tid_exons[line.rname];
        strand = exons[0].strand;
        # for exon in exons:
        #     print exon;

        cigar = line.SplitCigar()

        if (strand == '+'):
            cigar.reverse() #stack
            regInd = 0
            curr = line.pos
            #In which region does it start
            last = 0
            while last + regions[regInd][1] - regions[regInd][0] + 1 < curr:
                last += regions[regInd][1] - regions[regInd][0] + 1
                regInd = regInd + 1

            posOnRef = regions[regInd][0] + curr - last - 1

            newCigar = ""
            regionSize = regions[regInd][1] - regions[regInd][0] + 1

            # last = 0
            curr -= 1
            while cigar:
                c,op = cigar.pop()
                count = int(c)
                if op == 'I' or op == 'S' or op == 'H':
                    newCigar += str(count) + op
                    continue
                if count + curr > regionSize + last:
                    take = regionSize + last - curr
                    if regInd < len(regions) - 1:
                        if take > 0:
                            newCigar += str(take) + op
                        cigar.append((count - take, op))
                        newCigar += str(regions[regInd+1][0] - regions[regInd][1] - 1) + 'N'
                    else:
                        newCigar += str(c) + op
                        break

                    regInd = regInd + 1
                    last += regionSize
                    regionSize = regions[regInd][1] - regions[regInd][0] + 1
                    curr += take
                else:
                    curr += count
                    newCigar += str(count) + op

            while cigar:
                c,op = cigar.pop()
                newCigar += str(c) + op

            newLine = copy.deepcopy(line);
            newLine.cigar = newCigar
            newLine.rname = transToSeq[name][0]
            newLine.pos = posOnRef
            # if transToSeq[name][1] == '-':
            #     newLine.flag ^= 0x10
            newSam.append(newLine)

        elif (strand == '-'):
            #########################
            ### This is different ###
            trans_len = 0;
            for i in xrange(0, len(regions)):
                trans_len += regions[i][1] - regions[i][0] + 1;
            #########################

            ### Don't reverse the cigar.
            # cigar.reverse() #stack
            regInd = 0
            #########################
            ### This is different ###
            # curr = line.pos
            pos_end = line.pos + line.CalcReferenceLengthFromCigar() - 1;
            # print line.pos, line.CalcReferenceLengthFromCigar(), pos_end;
            curr = trans_len - pos_end;
            # print curr;
            #########################
            #In which region does it start
            last = 0
            while last + regions[regInd][1] - regions[regInd][0] + 1 < curr:
                last += regions[regInd][1] - regions[regInd][0] + 1
                regInd = regInd + 1

            posOnRef = regions[regInd][0] + curr - last;
            # print posOnRef;

#            newCigar = ""
            newCigarList = [];
            regionSize = regions[regInd][1] - regions[regInd][0] + 1

            # last = 0
            curr -= 1
            while cigar:
                c,op = cigar.pop()
                count = int(c)
                if op == 'I' or op == 'S' or op == 'H':
                    # newCigar += str(count) + op
                    newCigarList.append(str(count) + op);
                    continue
                if count + curr > regionSize + last:
                    take = regionSize + last - curr
                    if regInd < len(regions) - 1:
                        if take > 0:
                            # newCigar += str(take) + op
                            newCigarList.append(str(take) + op);
                        cigar.append((count - take, op))
#                        newCigar += str(regions[regInd+1][0] - regions[regInd][1] - 1) + 'N'
                        newCigarList.append(str(regions[regInd+1][0] - regions[regInd][1] - 1) + 'N');
                    else:
                        # newCigar += str(c) + op
                        newCigarList.append(str(c) + op);
                        break

                    regInd = regInd + 1
                    last += regionSize
                    regionSize = regions[regInd][1] - regions[regInd][0] + 1
                    curr += take
                else:
                    curr += count
                    # newCigar += str(count) + op
                    newCigarList.append(str(count) + op);

            while cigar:
                c,op = cigar.pop()
                # newCigar += str(c) + op
                newCigarList.append(str(c) + op);

            newLine = copy.deepcopy(line);
            # print 'Tu sam!';
#            newCigar = ''.join(newCigarList[::-1]);
            newCigar = ''.join(newCigarList[:]);
            newLine.cigar = newCigar;
            newLine.flag = newLine.flag ^ 0x10; # Change the reverse flag.
            newLine.rname = transToSeq[name][0]
            newLine.seq = ReverseComplement(newLine.seq);
            newLine.pos = posOnRef
            # if transToSeq[name][1] == '-':
            #     newLine.flag ^= 0x10
            newSam.append(newLine)
            
    return newSam

def ReverseComplement(seq):
    l = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'};
    return ''.join([l[s] for s in seq[::-1]]);

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "Error - 3 arguments needed: in.sam in.gtf out.sam"
        sys.exit(0)

    headers, sam_lines = utility_sam.LoadSAM(sys.argv[1])
    gtf = open(sys.argv[2])
    tid_exons, transToSeq = transExtract.parse(gtf)
    gtf.close()
    tid_regions = transExtract.makeRegions(tid_exons)

    newSam = solve(sam_lines, transToSeq, tid_regions, tid_exons)
    samOut = open(sys.argv[3], "w")
    for head in headers:
        samOut.write(head + '\n')
    for sam in newSam:
        samOut.write("\t".join([sam.qname, str(sam.flag), sam.rname, str(sam.pos), \
            str(sam.mapq), sam.cigar, sam.mrnm, str(sam.mpos), str(sam.isize), sam.seq, sam.qual]) + "\n");
    samOut.close();
