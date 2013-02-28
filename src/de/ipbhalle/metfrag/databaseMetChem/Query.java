/*
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/

package de.ipbhalle.metfrag.databaseMetChem;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.sql.Connection;
import java.sql.Driver;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.nonotify.NNMolecule;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import de.ipbhalle.metfrag.main.Config;
import de.ipbhalle.metfrag.massbankParser.Peak;
import de.ipbhalle.metfrag.spectrum.WrapperSpectrum;
import de.ipbhalle.metfrag.tools.PPMTool;

public class Query {
	
	private Connection con;
	private Map<String, Integer> databaseToID;
	String url;
    String username;
    String password;
    
    
    /**
	 * Instantiates a new query.
	 *
	 * @param username the username
	 * @param password the password
	 * @param url the url
	 */
	public Query(String username, String password, String url)
	{
		String driver = "org.postgresql.Driver"; 
		Driver driverPostgres = new org.postgresql.Driver();
		try {
			Class.forName(driver);
			DriverManager.registerDriver (driverPostgres);
	        //database data
	        this.url = url + "?protocolVersion=2";
	        this.username = username;
	        this.password = password;
	        con = DriverManager.getConnection(url, username, password);
	        
	        databaseToID = new HashMap<String, Integer>();
	        PreparedStatement pstmtDatabase = con.prepareStatement("SELECT library_id, library_name from library;");
	        ResultSet res = pstmtDatabase.executeQuery();
	        while(res.next()){
	        	databaseToID.put(res.getString(2), res.getInt(1));
	        }
	        
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	/*	finally
		{
			try {
				con.close();
				DriverManager.deregisterDriver(driverPostgres);
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}*/
		
	}
    
    
    /**
     * Gets the sorted candidates. This is the MassStruct method to retrieve candidates in
     * a sorted way. Reference: http://dx.doi.org/10.2390/biecoll-jib-2011-157
     *
     * @param ws the spectrum
     * @param searchPPM the exact mass restriction
     * @param mzabs the mzabs
     * @param mzppm the mzppm
     * @return the sorted candidates
     */
    public List<String> getSortedCandidates(WrapperSpectrum ws, String searchPPM, double mzabs, double mzppm) {
        List<String> SortedList = new ArrayList<String>();
        ResultSet rs = null;
        PreparedStatement pstmt = null;
        String query = null;
        try {
            //try to parse restriction in ppm
            double ppm = Double.parseDouble(searchPPM);
            query = "SELECT"
                    + " substance.accession, Score"
                    + " FROM substance, compound" + ", library,"
                    + " (SELECT inchi_key_1,"
                    + " COUNT(DISTINCT MCS.cluster_id) AS Score"
                    + " FROM substance,"
                    + " library,"
                    + " (SELECT MIN(compound_id) AS firstcompound_id"
                    + " FROM compound"
                    + " WHERE exact_mass BETWEEN " + calcMinValue(ws.getExactMass(), mzabs, mzppm)
                    + " AND " + calcMaxValue(ws.getExactMass(), mzabs, mzppm)
                    + " GROUP BY inchi_key_1) AS firstcandidate,"
                    + " compound AS Candidates"
                    + " LEFT OUTER JOIN (SELECT mcs.mcs_structure,"
                    + " mz_cluster.cluster_id"
                    + " FROM mcs,"
                    + " mz_cluster"
                    + " WHERE mcs.mz_cluster_id = mz_cluster.cluster_id"
                    + " AND (";
            for (Peak peak : ws.getPeakList()) {
                if (peak.equals(ws.getPeakList().lastElement())) {
                    query += "(min >= "
                            + calcMinValue(peak.getMass(), mzabs, mzppm)
                            + " AND max <= "
                            + calcMaxValue(peak.getMass(), mzabs, mzppm)
                            + "))) AS MCS "
                            + " ON (MCS.mcs_structure <= Candidates.mol_structure)"
                            + " WHERE substance.compound_id = firstcompound_id"
                            + " AND   substance.library_id = library.library_id"
                            + " AND   library_name = 'pubchem'"
                            + " AND   Candidates.compound_id = firstcompound_id"
                            + " GROUP BY accession, inchi_key_1"
                            + " ORDER BY Score DESC) AS results"
                            + " WHERE exact_mass BETWEEN " + calcMinValue(ws.getExactMass(), mzabs, mzppm)
                            + " AND " + calcMaxValue(ws.getExactMass(), mzabs, mzppm)
                            + " AND substance.compound_id = compound.compound_id"
                            + " AND substance.library_id = library.library_id"
                            + " AND library_name = 'pubchem'"
                            + " AND compound.inchi_key_1 = results.inchi_key_1;";
                } else {
                    query += "(min >= "
                            + calcMinValue(peak.getMass(), mzabs, mzppm)
                            + " AND max <= "
                            + calcMaxValue(peak.getMass(), mzabs, mzppm)
                            + ") OR ";
                }
            }
        } catch (Exception e) {
            //it is not a ppm but a sum formula
            query = "SELECT"
                    + " substance.accession, Score"
                    + " FROM substance, compound" + ", library,"
                    + " (SELECT inchi_key_1,"
                    + " COUNT(DISTINCT MCS.cluster_id) AS Score"
                    + " FROM substance,"
                    + " library,"
                    + " (SELECT MIN(compound_id) AS firstcompound_id"
                    + " FROM compound"
                    + " WHERE formula = " + searchPPM
                    + " GROUP BY inchi_key_1) AS firstcandidate,"
                    + " compound AS Candidates"
                    + " LEFT OUTER JOIN (SELECT mcs.mcs_structure,"
                    + " mz_cluster.cluster_id"
                    + " FROM mcs,"
                    + " mz_cluster"
                    + " WHERE mcs.mz_cluster_id = mz_cluster.cluster_id"
                    + " AND (";
            for (Peak peak : ws.getPeakList()) {
                if (peak.equals(ws.getPeakList().lastElement())) {
                    query += "(min >= "
                            + calcMinValue(peak.getMass(), mzabs, mzppm)
                            + " AND max <= "
                            + calcMaxValue(peak.getMass(), mzabs, mzppm)
                            + "))) AS MCS "
                            + " ON (MCS.mcs_structure <= Candidates.mol_structure)"
                            + " WHERE substance.compound_id = firstcompound_id"
                            + " AND   substance.library_id = library.library_id"
                            + " AND   library_name = 'pubchem'"
                            + " AND   Candidates.compound_id = firstcompound_id"
                            + " GROUP BY accession, inchi_key_1"
                            + " ORDER BY Score DESC) AS results"
                            + " WHERE formula = " + searchPPM
                            + " AND substance.compound_id = compound.compound_id"
                            + " AND substance.library_id = library.library_id"
                            + " AND library_name = 'pubchem'"
                            + " AND compound.inchi_key_1 = results.inchi_key_1;";
                } else {
                    query += "(min >= "
                            + calcMinValue(peak.getMass(), mzabs, mzppm)
                            + " AND max <= "
                            + calcMaxValue(peak.getMass(), mzabs, mzppm)
                            + ") OR ";
                }
            }
        }
        try {
        	con = DriverManager.getConnection(url, username, password);
            pstmt = con.prepareStatement(query);
            rs = pstmt.executeQuery();
        } catch (Exception e) {
            System.err.println("Problems with Statement: "+e.getLocalizedMessage());
        }
        try {
            while (rs.next()) {
                SortedList.add(rs.getString(1));
            }
        } catch (Exception e) {
            System.err.println("Problems with Resultset: " + e.getLocalizedMessage());
        }
        return SortedList;
    }


    private Double calcMinValue(Double d, double mz_abs, double ppm) {
        return d - PPMTool.getPPMDeviation(d, ppm) - mz_abs;
    }

    private Double calcMaxValue(Double d, double mz_abs, double ppm) {
        return d + PPMTool.getPPMDeviation(d, ppm) + mz_abs;
    } 
    
	
	/**
	 * 
	 * @param lowerBound
	 * @param upperBound
	 * @param database
	 * @param uniqueInchi
	 * @return
	 */
	public List<CandidateMetChem> queryByMass(double lowerBound, double upperBound, String database, boolean uniqueInchi)
	{
		if(this.databaseToID.get(database) == null)
		{
			String possible = "";
			for (String db : this.databaseToID.keySet()) {
				possible += "\"" + db + "\" ";
			}
			System.err.println("Wrong database name given! Possible are: " + possible);
			return null;
		}
		
		List<CandidateMetChem> candidatesString = new ArrayList<CandidateMetChem>(); 
		try
		{
			con = DriverManager.getConnection(url, username, password);
			con.setAutoCommit(false);
			
			//SELECT c.compound_id from compound c inner join substance s on s.compound_id = c.compound_id 
			//where exact_mass >= 272.071 and exact_mass <= 272.079 and s.library_id = 2;
			
		    PreparedStatement pstmt = null;
		    if(!uniqueInchi) pstmt = con.prepareStatement("SELECT c.compound_id, accession from compound c inner join substance s on s.compound_id = c.compound_id " +
		    		"where exact_mass between ? and ? and s.library_id = ?;");
		    else pstmt = con.prepareStatement("SELECT distinct on (c.inchi_key_1) c.compound_id, accession from compound c inner join substance s on s.compound_id = c.compound_id " +
		    		"where exact_mass between ? and ? and s.library_id = ?;");
		   
		    pstmt.setDouble(1, lowerBound);
		    pstmt.setDouble(2, upperBound);
		    pstmt.setInt(3, this.databaseToID.get(database));
		    
		    //System.out.println(pstmt.toString());
		    ResultSet res = pstmt.executeQuery();
	       
	        while(res.next()){
	        	candidatesString.add(new CandidateMetChem(res.getInt(1), res.getString(2)));
	        }
		}
		catch(SQLException e)
		{
			System.err.println("SQL error! " + e.getMessage());
			e.printStackTrace();
		}
/*		finally
		{
			try {
				con.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}*/
		
        
        return candidatesString;
	}
	
	/**
	 * 
	 * @param formula
	 * @param database
	 * @param uniqueInchi
	 * @return
	 * @throws SQLException
	 */
	public List<CandidateMetChem> queryByFormula(String formula, String database, boolean uniqueInchi) throws SQLException
	{
		if(this.databaseToID.get(database) == null)
		{
			String possible = "";
			for (String db : this.databaseToID.keySet()) {
				possible += "\"" + db + "\" ";
			}
			System.err.println("Wrong database name given! Possible are: " + possible);
			return null;
		}
		
		List<CandidateMetChem> candidatesString = new ArrayList<CandidateMetChem>(); 
		
		try
		{
			con = DriverManager.getConnection(url, username, password);
			PreparedStatement pstmt = null;
			if(!uniqueInchi) pstmt = con.prepareStatement("SELECT c.compound_id, accession from compound c inner join substance s on s.compound_id = c.compound_id " +
					"where formula = ? and s.library_id = ?;");
			else pstmt = con.prepareStatement("SELECT distinct on (c.inchi_key_1) c.compound_id, accession from compound c inner join substance s on s.compound_id = c.compound_id " +
					"where formula = ? and s.library_id = ?;");
		    pstmt.setString(1, formula);
		    pstmt.setInt(2, this.databaseToID.get(database));
		    //System.out.println(pstmt.toString());
	        ResultSet res = pstmt.executeQuery();
	        while(res.next()){
	        	candidatesString.add(new CandidateMetChem(res.getInt(1), res.getString(2)));
	        }
	        pstmt.close();
		}
		catch(SQLException e)
		{
			System.err.println("SQL error! " + e.getMessage());
			e.printStackTrace();
		}
 /*       finally
		{
			try {
				con.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}*/
        
        return candidatesString;
	}
	
	
	/**
	 * Open connection.
	 */
	public void openConnection()
	{
		try {
			con = DriverManager.getConnection(url, username, password);
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public void closeConnection()
	{
		try {
			con.close();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Gets the compound with already opened connection
	 *
	 * @param compoundID the compound id
	 * @return the compound
	 * @throws SQLException the sQL exception
	 * @throws CDKException the cDK exception
	 */
	public IAtomContainer getCompoundConnectionOpened(Integer compoundID) throws SQLException, CDKException
	{
		IAtomContainer ret = null;
		
		PreparedStatement pstmt = con.prepareStatement("SELECT v3000(mol_structure) from compound c where compound_id = ? limit 1;");
	    pstmt.setInt(1, compoundID);
	    
        ResultSet res = pstmt.executeQuery();
        while(res.next()){
        	String molString = res.getString(1);
        	MDLV3000Reader reader = new MDLV3000Reader(new ByteArrayInputStream(molString.getBytes()));
        	ret = (IAtomContainer)reader.read(new NNMolecule());
        }
        pstmt.close();
        		
        return ret;
	}
	
	/**
	 * Gets the compound using identifier. (Previously opened connection)
	 *
	 * @param identifier the identifier
	 * @param database the database (pubchem, kegg, chebi, knapsack)
	 * @return the compound using identifier
	 * @throws Exception 
	 */
	public IAtomContainer getCompoundUsingIdentifierConnectionOpen(String identifier, String database) throws Exception
	{
		Integer databaseID = databaseToID.get(database);
		IAtomContainer ret = null;
		
		if(databaseID == null)
			throw new Exception("Database not supported!");
		
		try
		{
			PreparedStatement pstmt = con.prepareStatement("SELECT v3000(mol_structure) from compound c inner join substance s on s.compound_id = c.compound_id " +
			"where s.library_id = ? and s.accession = ?;");
		    pstmt.setInt(1, databaseID);
		    pstmt.setString(2, identifier);
		    //System.out.println(pstmt.toString());
	        ResultSet res = pstmt.executeQuery();
	        while(res.next()){
	        	String molString = res.getString(1);
	        	MDLV3000Reader reader = new MDLV3000Reader(new ByteArrayInputStream(molString.getBytes()));
	        	ret = (IAtomContainer)reader.read(new NNMolecule());
	        }
	        pstmt.close();
		}
		catch(SQLException e)
		{
			System.err.println("SQL error! " + e.getMessage());
			e.printStackTrace();
		}
        
        return ret;
	}
	
	
	
	
	/**
	 * Gets the compound.
	 *
	 * @param compoundID the compound id
	 * @return the compound
	 * @throws SQLException 
	 * @throws CDKException 
	 */
	public IAtomContainer getCompound(Integer compoundID) throws SQLException, CDKException
	{
		IAtomContainer ret = null;
		
		con = DriverManager.getConnection(url, username, password);
		PreparedStatement pstmt = con.prepareStatement("SELECT v3000(mol_structure) from compound c where compound_id = ? limit 1;");
	    pstmt.setInt(1, compoundID);
	    
        ResultSet res = pstmt.executeQuery();
        while(res.next()){
        	String molString = res.getString(1);
        	MDLV3000Reader reader = new MDLV3000Reader(new ByteArrayInputStream(molString.getBytes()));
        	ret = (IAtomContainer)reader.read(new NNMolecule());
        }
        pstmt.close();
        con.close();
		
        return ret;
	}
	
	
	/**
	 * Gets the compound using identifier.
	 *
	 * @param identifier the identifier
	 * @param database the database (pubchem, kegg, chebi, knapsack)
	 * @return the compound using identifier
	 * @throws Exception 
	 */
	public IAtomContainer getCompoundUsingIdentifier(String identifier, String database) throws Exception
	{
		Integer databaseID = databaseToID.get(database);
		IAtomContainer ret = null;
		
		if(databaseID == null)
			throw new Exception("Database not supported!");
		
		try
		{
			con = DriverManager.getConnection(url, username, password);
			PreparedStatement pstmt = con.prepareStatement("SELECT v3000(mol_structure) from compound c inner join substance s on s.compound_id = c.compound_id " +
			"where s.library_id = ? and s.accession = ?;");
		    pstmt.setInt(1, databaseID);
		    pstmt.setString(2, identifier);
		    //System.out.println(pstmt.toString());
	        ResultSet res = pstmt.executeQuery();
	        while(res.next()){
	        	String molString = res.getString(1);
	        	MDLV3000Reader reader = new MDLV3000Reader(new ByteArrayInputStream(molString.getBytes()));
	        	ret = (IAtomContainer)reader.read(new NNMolecule());
	        }
	        pstmt.close();
		}
		catch(SQLException e)
		{
			System.err.println("SQL error! " + e.getMessage());
			e.printStackTrace();
		}
  /*      finally
		{
			try {
				con.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}*/
        
        return ret;
	}
	
}
