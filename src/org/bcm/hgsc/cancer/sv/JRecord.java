package org.bcm.hgsc.cancer.sv;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.UUID;

import org.bcm.hgsc.cancer.utils.Orientation;

public class JRecord implements Serializable{

	/**
	 * serial version as of 14 March 2014
	 */
	private static final long serialVersionUID = 1L;
	private final String chra;
	private final String chrb;
	private final long posa;
	private final long posb;
	private final long posal;
	private final long posah;
	private final long posbl;
	private final long posbh;
	private final Orientation oapos;
	private final Orientation obpos;
	private final int ebases;
	private final String calleruuid;
	private final String uuid;
	//	private final boolean flipped;

	public void printJunction(BufferedWriter stream){
		boolean isDel = ebases < AppendSVGraph.minSVSize;
		try {
			stream.write(chra + "\t" + posa + "\t" + chrb + "\t" + posb + "\t" + oapos + "\t" + obpos + "\t" + (isDel ? 0 + "\tDEL" : ebases + "\tMIS") + "\t" + SVGraph.callers.get(this.calleruuid).getName() + "\n");
			stream.flush();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void printJunction(String fileName){
		try{
			if (fileName != null){
				printJunction(new BufferedWriter(new PrintWriter(new FileWriter(fileName, true))));
			} else {
				printJunction(new BufferedWriter(new OutputStreamWriter(System.out)));
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public JRecord(String chra, long posa, String chrb, long posb, Orientation oa, Orientation ob, int ebases, Caller caller){
		this.chra = chra;
		this.chrb = chrb;
		this.posa = posa;
		this.posb = posb;
		this.oapos = oa;
		this.obpos = ob;
		this.ebases = ebases;
		this.calleruuid = caller.getUUID();
		this.posal = posa - caller.getBuffer();
		this.posah = posa + caller.getBuffer();
		this.posbl = posb - caller.getBuffer();
		this.posbh = posb + caller.getBuffer();
		this.uuid = UUID.randomUUID().toString();
	}

	public JRecord(String chra, long posa, String chrb, long posb, Orientation oa, Orientation ob, int ebases, Caller caller, String uuid){
		this.chra = chra;
		this.chrb = chrb;
		this.posa = posa;
		this.posb = posb;
		this.oapos = oa;
		this.obpos = ob;
		this.ebases = ebases;
		this.calleruuid = caller.getUUID();
		this.posal = posa - caller.getBuffer();
		this.posah = posa + caller.getBuffer();
		this.posbl = posb - caller.getBuffer();
		this.posbh = posb + caller.getBuffer();
		this.uuid = uuid;
	}

	public JRecord(JRecord j){
		this.chra = j.chra;
		this.chrb = j.chrb;
		this.posa = j.posa;
		this.posb = j.posb;
		this.oapos = j.oapos;
		this.obpos = j.obpos;
		this.ebases = j.ebases;
		this.calleruuid = j.calleruuid;
		this.posal = j.posal;
		this.posah = j.posah;
		this.posbl = j.posbl;
		this.posbh = j.posbh;
		this.uuid = j.uuid;
	}

	public String getUUID(){
		return this.uuid;
	}

	public String chra() {
		return this.chra;
	}

	public long posa() {
		return this.posa;
	}

	public Orientation oapos() {
		return this.oapos;
	}

	public String chrb() {
		return this.chrb;
	}

	public long posb() {
		return this.posb;
	}

	public Orientation obpos() {
		return this.obpos;
	}

	public int ebases() {
		return this.ebases;
	}

	public long getPosal() {
		return posal;
	}

	public long getPosah() {
		return posah;
	}

	public long getPosbl() {
		return posbl;
	}

	public long getPosbh() {
		return posbh;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((calleruuid == null) ? 0 : calleruuid.hashCode());
		result = prime * result + ((chra == null) ? 0 : chra.hashCode());
		result = prime * result + ((chrb == null) ? 0 : chrb.hashCode());
		result = prime * result + ebases;
		result = prime * result + ((oapos == null) ? 0 : oapos.hashCode());
		result = prime * result + ((obpos == null) ? 0 : obpos.hashCode());
		result = prime * result + (int) (posa ^ (posa >>> 32));
		result = prime * result + (int) (posb ^ (posb >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof JRecord)) {
			return false;
		}
		JRecord other = (JRecord) obj;
		if (calleruuid == null) {
			if (other.calleruuid != null) {
				return false;
			}
		} else if (!calleruuid.equals(other.calleruuid)) {
			return false;
		}
		if (chra == null) {
			if (other.chra != null) {
				return false;
			}
		} else if (!chra.equals(other.chra)) {
			return false;
		}
		if (chrb == null) {
			if (other.chrb != null) {
				return false;
			}
		} else if (!chrb.equals(other.chrb)) {
			return false;
		}
		if (ebases != other.ebases) {
			return false;
		}
		if (oapos != other.oapos) {
			return false;
		}
		if (obpos != other.obpos) {
			return false;
		}
		if (posa != other.posa) {
			return false;
		}
		if (posb != other.posb) {
			return false;
		}
		return true;
	}

	public String getCallerUUID(){
		return this.calleruuid;
	}


}
