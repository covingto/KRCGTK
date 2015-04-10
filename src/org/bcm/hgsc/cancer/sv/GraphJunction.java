package org.bcm.hgsc.cancer.sv;

import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;

import org.bcm.hgsc.cancer.utils.Orientation;

public class GraphJunction extends JRecord{
	public static enum HitType {
		Full, Left, Right, None
	}
	/**
	 * serial version as of 14 March 2014
	 */
	private static final long serialVersionUID = 1L;
	private int familyInt = 0;
	private int familyCount = 0;
	public final Set<String> partners = new HashSet<String>();
	public final Set<String> partnerCallers = new HashSet<String>();
	// =====================
	public static final Comparator<GraphJunction> aComparator = new Comparator<GraphJunction>(){

		@Override
		public int compare(GraphJunction o1, GraphJunction o2) {
			if (o1.chra().equals(o2.chra())){
				return (int) (o1.posa() - o2.posa());
			} else {
				return o1.chra().compareTo(o2.chra());
			}
		}
		
	};
	
	// =====================
	public static final Comparator<GraphJunction> bComparator = new Comparator<GraphJunction>(){

		@Override
		public int compare(GraphJunction o1, GraphJunction o2) {
			if (o1.chrb().equals(o2.chrb())){
				return (int) (o1.posb() - o2.posb());
			} else {
				return o1.chrb().compareTo(o2.chrb());
			}
		}
		
	};
	
	// =====================
	public static final Comparator<GraphJunction> interactionCountComparator = new Comparator<GraphJunction>(){
		@Override
		public int compare(GraphJunction o1, GraphJunction o2){
			return o1.partnerCallers.size() - o2.partnerCallers.size();
			
		}
	};
	
	// =====================
	public GraphJunction(String chra, long posa, String chrb, long posb,
			Orientation oa, Orientation ob, int ebases, Caller caller) {
		super(chra, posa, chrb, posb, oa, ob, ebases, caller);
	}
	
	public GraphJunction(String chra, long posa, String chrb, long posb,
			Orientation oa, Orientation ob, int ebases, Caller caller, String uuid) {
		super(chra, posa, chrb, posb, oa, ob, ebases, caller, uuid);
	}
	
	public GraphJunction(JRecord j){
		super(j);
	}

	public void addPartner(String p){
		partners.add(p);
	}
	
	public void addPartnerCaller(String c){
		partnerCallers.add(c);
	}

	/**
	 * Returns a comparison of the junctions.  If not a match return the comparason of a to a.
	 */
	public boolean intersects(GraphJunction arg0) {
		if (this.getCallerUUID().equals(arg0.getCallerUUID())){ return false; } // while this may actually intersect in base space, callers can't validate themselves so we return false
		boolean compab = compareJa(arg0.chrb(), arg0.getPosbl(), arg0.getPosbh());
		if (compab){ return true; }
		boolean compba = compareJb(arg0.chra(), arg0.getPosal(), arg0.getPosah());
		if (compba){ return true; }
		boolean compbb = compareJb(arg0.chrb(), arg0.getPosbl(), arg0.getPosbh());
		if (compbb){ return true; }
		return compareJa(arg0.chra(), arg0.getPosal(), arg0.getPosah());
		
	}
	
	// ============================
	public HitType hitType(GraphJunction arg0){
		if (this.getCallerUUID().equals(arg0.getCallerUUID())){ return HitType.None; } // while this may actually intersect in base space, callers can't validate themselves so we return false
		final boolean compab = compareJa(arg0.chrb(), arg0.getPosbl(), arg0.getPosbh());
		final boolean compba = compareJb(arg0.chra(), arg0.getPosal(), arg0.getPosah());
		final boolean compbb = compareJb(arg0.chrb(), arg0.getPosbl(), arg0.getPosbh());
		final boolean compaa = compareJa(arg0.chra(), arg0.getPosal(), arg0.getPosah());
		if ((compab && compba) || (compba && compab) || (compaa && compbb)){
			return HitType.Full;
		} else if (compab || compaa ){ 
			return HitType.Left;
		} else if (compbb || compba){
			return HitType.Right;
		} else {
			return HitType.None;
		}
	}
	
	// ============================
	private boolean compareJa(String chr, long low, long high){
		if (this.chra().equals(chr)){
			if (this.getPosah() > low && this.getPosal() < high){ return true; }
			else { return false; }
		} else {
			return this.chra().equals(chr);
		}
	}
	
	private boolean compareJb(String chr, long low, long high){
		if (this.chrb().equals(chr)){
			if (this.getPosbh() > low && this.getPosbl() < high){ return true; }
			else { return false; }
		} else {
			return this.chrb().equals(chr);
		}
	}

	public int partnerCount() {
		return partners.size();
	}

	public void setFamilyInt(int familyInt) {
		this.familyInt = familyInt;
	}
	
	public void setFamilyCount(int familyCount){
		this.familyCount = familyCount;
	}
	
	public int getFamilyInt(){
		return this.familyInt;
	}
	
	public int getFamilyCount(){
		return this.familyCount;
	}
	
}
