package testcases;

import static org.junit.Assert.*;

import org.bcm.hgsc.utils.AlleleResolver.AlleleSet;
import org.junit.Test;

public class TestAlleleResolver {

	@Test
	public void testGetRightOffset() {
		assertEquals("Failed at simple SNP", 2, 
				AlleleSet.getRightOffset(
						new String("ATCG").getBytes(),
						new String("AACG").getBytes()));
		assertEquals("Failed at simple INDEL", 2,
				AlleleSet.getRightOffset(new String("ATCG").getBytes(), new String("ACG").getBytes()));
		assertEquals("Failed at repeat INDEL", 4,
				AlleleSet.getRightOffset(new String("ATCGCGAT").getBytes(), new String("ATCGAT").getBytes()));
		assertEquals("Failed at complex INDEL", 2,
				AlleleSet.getRightOffset(
						new String("ATCGCGAT").getBytes(),
						new String("ATAAT").getBytes()));
	}
	
	@Test
	public void testCorrectLeftOffsetRepeat(){
		byte[] a = new String("ATCGCGCGAT").getBytes();
		byte[] b = new String("ATCGAT").getBytes();
		int left = AlleleSet.getLeftOffset(a, b);
		int right = AlleleSet.getRightOffset(a, b);
		assertEquals("Left offset is not correct", 4, left);
		assertEquals("Right offset is not correct", 4, right);
		assertEquals("Correction not correct", 2, AlleleSet.correctLeftOffset(a, b, left, right));
	}
	
	@Test
	public void testCorrectLeftOffsetSimple(){
		byte[] a = new String("ATGGCG").getBytes();
		byte[] b = new String("ACG").getBytes();
		int left = AlleleSet.getLeftOffset(a, b);
		int right = AlleleSet.getRightOffset(a, b);
		assertEquals("Left offset is not correct", 1, left);
		assertEquals("Right offset is not correct", 2, right);
		assertEquals("Correction not correct", 1, AlleleSet.correctLeftOffset(a, b, left, right));
	}

	@Test
	public void testGetLeftOffset() {
		assertEquals("Failed at simple SNP", 1, 
				AlleleSet.getLeftOffset(
						new String("ATCG").getBytes(),
						new String("AACG").getBytes()));
		assertEquals("Failed at simple INDEL", 1,
				AlleleSet.getLeftOffset(new String("ATCG").getBytes(), new String("ACG").getBytes()));
		assertEquals("Failed at repeat INDEL", 4,
				AlleleSet.getLeftOffset(new String("ATCGCGAT").getBytes(), new String("ATCGAT").getBytes()));
		assertEquals("Failed at complex INDEL", 2,
				AlleleSet.getLeftOffset(
						new String("ATCGCGAT").getBytes(),
						new String("ATAAT").getBytes()));
	}
	
	@Test
	public void testSimplifyAllele_simple(){
		byte[] a = new String("ATCG").getBytes();
		byte[] b = new String("AACG").getBytes();
		int left = AlleleSet.getLeftOffset(a, b);
		int right = AlleleSet.getRightOffset(a, b);
		int left_c = AlleleSet.correctLeftOffset(a, b, left, right);
		assertEquals("Failed to get left offset", 1, left);
		assertEquals("Failed to get right offset", 2, right);
		assertEquals("Failed to get corrected left", 1, left_c);
		
	}
	
	@Test
	public void testSimplifyAlleleProblemTestCases(){
		// this generated the wrong output at some point
		// returned TC and TCCCCT when it should have returned C and CCCCT
		byte[] a = new String("TC").getBytes();
		byte[] b = new String("TCCCCT").getBytes();
		int left = AlleleSet.getLeftOffset(a, b);
		int right = AlleleSet.getRightOffset(a, b);
		int left_c = AlleleSet.correctLeftOffset(a, b, left, right);
		assertEquals("Failed to get left offset", 1, left);
		assertEquals("Failed to get right offset", 0, right);
		assertEquals("Failed to get corrected left", 1, left_c);
		// solution, the default for @ref AlleleSet.getLeftOffset() was sto return 
		// 0 if there were no matches.  This should be to return the index position of the last
		// available base (minLength - 1) as the offset.  Fixed and re-tested.  All tests pass 
		// KRC 10Nov2015
	}

}
