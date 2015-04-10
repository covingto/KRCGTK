package org.bcm.hgsc.cancer.pacbio;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class RefMaker {
	
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		if (args.length < 3) {
			// no arguments indicated
			System.out.println("RefMaker.jar [reffasta] [blastout] [outputdirectory]");
			System.exit(10);
		}
		System.out.println("Starting analysis");
		File reffasta = new File(args[0]);
		File splitreads = new File(args[1]);
		String outbase = args[2];
		// File logoutput = new File(outbase + ".log");
		File assemblies = new File(outbase, "assemblies");
		assemblies.mkdirs(); // make the directories if they don't exist
		
		IndexedFastaSequenceFile fasta = new IndexedFastaSequenceFile(reffasta);
		
		System.out.println("Building read library");
		BlastStruct reads = new BlastStruct(splitreads, fasta);
		
		System.out.println("Filtering reads");
		BlastGroup.filterUnmapped(reads, new File(outbase + ".unmapped"));
		
		System.out.println("Building and writing assemblies");
		makeAssemblies(reads, assemblies, reffasta);

	}
	
	
	
	public static void makeAssemblies(BlastStruct reads, File assembliesDir, File reffasta) throws IOException {
		// here we need to make the assemblies for the reads
		for (String k : reads.getMappings().keySet()){
			System.out.println("Processing read " + k);
			// k is the key for this set
			BlastGroup bg = reads.getMappings().get(k);
			List<SyntheticReference> srefs = new ArrayList<SyntheticReference>();
			for (int i = 0; i < bg.numMappings(); i++){
				buildRef(new BlankSyntheticReference(), i, bg, 0, 0, srefs);
			}
			if (srefs.size() < 1) { continue; } // don't write out empty records
			// write the references to disk
			FileWriter fw = new FileWriter(new File(assembliesDir, "synthetic_" + k.replace("/", "_") + ".fasta"));
			BufferedWriter writer = new BufferedWriter(fw);
			for (SyntheticReference sref : srefs){
				writer.write(sref.fastaName());
				writer.newLine();
				writer.write(sref.fastaSeq());
				writer.newLine();
			}
			fw.close();
		}
	}
	
	private static void buildRef(
			SyntheticReference reforg, // the reference that we are building onto
			Integer index, // the index of the next blast row
			BlastGroup bg, // the blast group object
			Integer lastposq, // the last position of the query covered
			Integer lastposr, // the last position of the growing reference (nqend - start)
			List<SyntheticReference> srefs // the list of references that we build onto with new results
			) {
		// we build onto reforg, which is the synthetic reference up to this point
		BlastRow br = bg.getMapping(index); // this is the next row to get
		if (br.unique == BlastRow.REDUNDANT){ return; } // we don't add redundant values
		System.out.println(br.toString());
		Integer nqstart = br.qstart;
		Integer nqend = br.qend;
		Integer length = bg.length();
		
		// set the clips, clips are from the perspective of the next read to be incorporated
		// leftclip is left of the br and right clip is right of br with respect to the query read
		
		Integer leftclip = lastposq == 0 ? nqstart : -(( (lastposq + nqstart) / 2 ) - nqstart);
		Integer rightclip = length - nqend;
		System.out.println("Move to the left: " + leftclip);
		System.out.println("Move to the right: " + rightclip);
		// clip the growing reference (reforg) and the next sequence
		SyntheticReference leftsr;
		SyntheticReference newsr;
		try{
			leftsr = reforg.clip(leftclip);
			// make the next sequence
			newsr = leftsr.add(br, leftclip, rightclip);
			
		} catch (ExcessiveClipException e) {
			System.out.println(e.getMessage());
			return;
		} catch (DoubleClipException e) {
			// TODO Auto-generated catch block
			System.out.println(e.getMessage());
			return;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
		
		
		
		/*
		Integer cliplen = lastposq == 0 ? -nqstart : (( (lastposq + nqstart) / 2 ) - nqstart);
		Integer leftclip = lastposr == 0 ? 0 : lastposr - cliplen;
		Integer rstart = br.sstart + (cliplen * o);
		Integer tail = (bg.length()-br.qend);
		Integer rend = br.send + (tail * o);
		System.out.println("LastposQuery=" + lastposq + ", LastposSR=" + lastposr + ", NQStart=" + nqstart + ", Cliplen=" + cliplen + ", NQEnd=" + nqend + ", Length=" + bg.length() + 
				" NSStart=" + br.sstart + ", NSEnd=" + br.send);
		
		if (lastposq == 0){
			// this is the first read in this batch, so we don't have a current reference, we just use the blank one
			leftQuery = reforg;
		} else {
			// we have a read already and need to do a chewback to the leftoffset position
			leftQuery = reforg.chewBack(leftclip);
		}
		
		rightQuery = new SyntheticReference(br.schr + "_" + rstart + "_" + rend + "_" + br.o, 
				getRef(reffasta, br.schr, rstart, rend, br.o));
		SyntheticReference newref = SyntheticReference.add(leftQuery, rightQuery);
		*/
		
		srefs.add(newsr);
		
		for (int i = index + 1; i < bg.numMappings(); i++){
			/*
			private static void buildRef(
					SyntheticReference reforg, // the reference that we are building onto
					Integer index, // the index of the next blast row
					BlastGroup bg, // the blast group object
					Integer lastposq, // the last position of the query covered
					Integer lastposr, // the last position of the growing reference (nqend - start)
					IndexedFastaSequenceFile reffasta, // the fasta for the reference
					List<SyntheticReference> srefs // the list of references that we build onto with new results
					)
			*/
			buildRef(newsr, i, bg, nqend, newsr.fastaSeq().length() - rightclip, srefs);
		}
	}
	
	/*
	public static String getRef(
			IndexedFastaSequenceFile reforg, 
			String contig, 
			Integer refseqstart, 
			Integer refseqend, 
			BlastRow.Orientation o
			){
		System.out.println("Query reference " + contig + " " + refseqstart + " " + refseqend + " " + o);
		switch (o) {
			case POS:
				return new String(reforg.getSubsequenceAt(contig, refseqstart, refseqend).getBases());
			case NEG:
				return SequenceUtil.reverseComplement(new String(reforg.getSubsequenceAt(contig, refseqend, refseqstart).getBases()));
			default:
				System.out.println("Orientation for read not in positive or negative orientation!");
				return new String(reforg.getSubsequenceAt(contig, refseqstart, refseqend).getBases());
		}
	}
	*/
	
	/*
	public static String chewBack(
			IndexedFastaSequenceFile reffasta, 
			String contig, 
			Integer refstart, 
			Integer refend, 
			BlastRow.Orientation o,
			Integer qstart, 
			Integer qend, 
			Integer leftoffset,
			Integer rightoffset
			){
		Integer refseqstart = null;
		Integer refseqend = null;
		switch (o) {
			case POS:
				refseqstart = refstart - leftoffset;
				refseqend = refend + rightoffset;
				return reffasta.getSubsequenceAt(contig, refseqstart, refseqend).toString();
			case NEG:
				refseqstart = refend - leftoffset;
				refseqend = refstart + rightoffset;
				ReferenceSequence rs = reffasta.getSubsequenceAt(contig, refseqstart, refseqend);
				return SequenceUtil.reverseComplement(rs.toString());
			default:
				System.out.println("Orientation for read not in positive or negative orientation!");
				return reffasta.getSubsequenceAt(contig, refstart, refend).toString();
		}
	}
	*/

	

}
