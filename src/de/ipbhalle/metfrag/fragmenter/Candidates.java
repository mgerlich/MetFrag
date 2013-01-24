/*
*
* Copyright (C) 2009-2010 IPB Halle, Sebastian Wolf
*
* Contact: swolf@ipb-halle.de
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
*
*/


package de.ipbhalle.metfrag.fragmenter;

import java.io.StringReader;
import java.rmi.RemoteException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import net.sf.jniinchi.INCHI_RET;
import net.sf.jniinchi.JniInchiException;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import de.ipbhalle.metfrag.chemspiderClient.ChemSpider;
import de.ipbhalle.metfrag.keggWebservice.KeggRestService;
import de.ipbhalle.metfrag.keggWebservice.KeggWebservice;
import de.ipbhalle.metfrag.molDatabase.KEGGLocal;
import de.ipbhalle.metfrag.molDatabase.PubChemLocal;
import de.ipbhalle.metfrag.pubchem.PubChemWebService;
import de.ipbhalle.metfrag.tools.PPMTool;

public class Candidates {
	
	
	/**
	 * Gets the candidates online using the webservice interface from the databases
	 * 
	 * @param database the database
	 * @param databaseID the database id
	 * @param molecularFormula the molecular formula
	 * @param exactMass the exact mass
	 * @param searchPPM the search ppm
	 * @param useIPBProxy the use ipb proxy
	 * 
	 * @return the online
	 * 
	 * @throws Exception the exception
	 */
	public static Vector<String> getOnline(String database, String databaseID, String molecularFormula, double exactMass, double searchPPM, boolean useIPBProxy, PubChemWebService pubchem) throws Exception
	{
		Vector<String> candidates = new Vector<String>();
		
		if(database.equals("kegg") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!molecularFormula.isEmpty())
				//candidates = KeggWebservice.KEGGbySumFormula(molecularFormula);
				candidates = KeggRestService.KEGGbySumFormula(molecularFormula);
			else
				//candidates = KeggWebservice.KEGGbyMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
				candidates = KeggRestService.KEGGbyMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
		}
		else if(database.equals("chemspider") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!molecularFormula.isEmpty())
				candidates = ChemSpider.getChemspiderBySumFormula(molecularFormula);
			else
				candidates = ChemSpider.getChemspiderByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
		}
		else if(database.equals("pubchem") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!molecularFormula.isEmpty())
				candidates = pubchem.getHitsbySumFormula(molecularFormula, useIPBProxy);
			else
				candidates = pubchem.getHitsByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)), Integer.MAX_VALUE, useIPBProxy);
		}
		else if (!databaseID.equals(""))
		{	
			candidates = new Vector<String>();
			candidates.add(databaseID);
//			String[] idList = databaseID.split(",");
//			for (int i = 0; i < idList.length; i++) {
//				candidates.add(idList[i].trim());
//			}
		}
		
		return candidates;
	}
	
	/**
	 * Gets the candidates online using the webservice interface from the databases
	 * 
	 * @param database the database
	 * @param databaseID the database id
	 * @param molecularFormula the molecular formula
	 * @param exactMass the exact mass
	 * @param searchPPM the search ppm
	 * @param useIPBProxy the use ipb proxy
	 * 
	 * @return the online
	 * 
	 * @throws Exception the exception
	 */
	public static Vector<String> getOnline(String database, String databaseID, String molecularFormula, double exactMass, double searchPPM, 
			boolean useIPBProxy, PubChemWebService pubchem, boolean uniqueInchi, String chemspiderToken) throws Exception
	{
		Vector<String> candidates = new Vector<String>();
		
		if(database.equals("kegg") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!molecularFormula.isEmpty())
				//candidates = KeggWebservice.KEGGbySumFormula(molecularFormula);
				candidates = KeggRestService.KEGGbySumFormula(molecularFormula);
			else
				//candidates = KeggWebservice.KEGGbyMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
				candidates = KeggRestService.KEGGbyMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
			if(uniqueInchi) 
				candidates = removeDuplicatesByInchi(candidates, database, chemspiderToken);
		}
		else if(database.equals("chemspider") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!uniqueInchi) {
				if(!molecularFormula.isEmpty())
					candidates = ChemSpider.getChemspiderBySumFormula(molecularFormula);
				else
					candidates = ChemSpider.getChemspiderByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
			}
			else {
				if(!molecularFormula.isEmpty())
					candidates = ChemSpider.getChemspiderBySumFormulaUniqueInchi(molecularFormula, chemspiderToken);
				else
					candidates = ChemSpider.getChemspiderByMassUniqueInchi(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)), chemspiderToken);
			}
		}
		else if(database.equals("pubchem") && databaseID.equals(""))
		{
			//if(molecularFormula != "")
			if(!molecularFormula.isEmpty())
				candidates = pubchem.getHitsbySumFormula(molecularFormula, useIPBProxy, uniqueInchi);
			else
				candidates = pubchem.getHitsByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)), Integer.MAX_VALUE, useIPBProxy, uniqueInchi);
		}
		else if (!databaseID.equals(""))
		{	
			candidates = new Vector<String>();
			candidates.add(databaseID);
		}
		
		return candidates;
	}
	
	/**
	 * Gets the candidates locally using a local database.
	 * 
	 * @param database the database
	 * @param exactMass the exact mass
	 * @param searchPPM the search ppm
	 * @param databaseUrl the database url
	 * @param username the username
	 * @param password the password
	 * 
	 * @return the locally
	 * 
	 * @throws SQLException the SQL exception
	 * @throws ClassNotFoundException the class not found exception
	 * @throws RemoteException the remote exception
	 */
	public static List<String> getLocally(String database, double exactMass, double searchPPM, String databaseUrl, String username, String password) throws SQLException, ClassNotFoundException, RemoteException
	{
		List<String> candidates = new ArrayList<String>();
		
		if(database.equals("kegg"))
		{
			KEGGLocal kl = new KEGGLocal(databaseUrl, username, password);
			double deviation = PPMTool.getPPMDeviation(exactMass, searchPPM);
			candidates = kl.getHits(Integer.MAX_VALUE, (exactMass - deviation) , (exactMass + deviation));				
		}
		else if(database.equals("chemspider"))
		{
			candidates = ChemSpider.getChemspiderByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
		}
		else if(database.equals("pubchem"))
		{
			PubChemLocal pl = new PubChemLocal(databaseUrl, username, password);
			double deviation =PPMTool.getPPMDeviation(exactMass, searchPPM);
			candidates = pl.getHits((exactMass - deviation), (exactMass + deviation));
		}
		
		return candidates;
	}


	/**
	 * Gets the candidates locally using a local database.
	 * 
	 * @param database the database
	 * @param exactMass the exact mass
	 * @param searchPPM the search ppm
	 * @param databaseUrl the database url
	 * @param username the username
	 * @param password the password
	 * 
	 * @return the locally
	 * 
	 * @throws SQLException the SQL exception
	 * @throws ClassNotFoundException the class not found exception
	 * @throws RemoteException the remote exception
	 * @throws CDKException 
	 * @throws JniInchiException 
	 */
	public static List<String> getLocally(String database, double exactMass, double searchPPM, String databaseUrl, String username, 
			String password, boolean uniqueInchi, String chemspiderToken) throws SQLException, ClassNotFoundException, RemoteException, CDKException, JniInchiException
	{
		List<String> candidates = new ArrayList<String>();
		
		if(database.equals("kegg"))
		{
			KEGGLocal kl = new KEGGLocal(databaseUrl, username, password);
			double deviation = PPMTool.getPPMDeviation(exactMass, searchPPM);
			candidates = kl.getHits(Integer.MAX_VALUE, (exactMass - deviation) , (exactMass + deviation));		
			if(uniqueInchi) 
				candidates = removeDuplicatesByInchiLocally(candidates, database, databaseUrl, username, password, chemspiderToken);		
		}
		else if(database.equals("chemspider"))
		{
			if(!uniqueInchi) {
				candidates = ChemSpider.getChemspiderByMass(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)));
			}
			else {
				candidates = ChemSpider.getChemspiderByMassUniqueInchi(exactMass, (PPMTool.getPPMDeviation(exactMass, searchPPM)), chemspiderToken);
			}
		}
		else if(database.equals("pubchem"))
		{
			PubChemLocal pl = new PubChemLocal(databaseUrl, username, password);
			double deviation =PPMTool.getPPMDeviation(exactMass, searchPPM);
			candidates = pl.getHits((exactMass - deviation), (exactMass + deviation), uniqueInchi);
		}
		
		return candidates;
	}

	
	/**
	 * Gets the compound using the webservice.
	 * 
	 * @param database the database
	 * @param candidate the candidate
	 * @param pw the pw
	 * 
	 * @return the compound
	 * @throws RemoteException 
	 * @throws CDKException 
	 */
	public static IAtomContainer getCompound(String database, String candidate, 
			PubChemWebService pw, String chemspiderToken) throws RemoteException, CDKException
	{
		IAtomContainer molecule = null;
		
		if(database.equals("kegg"))
		{	
			if(candidate.startsWith("cpd:"))
				candidate = candidate.substring(4);
			
			//String candidateMol = KeggWebservice.KEGGgetMol(candidate, "");
			String candidateMol = KeggRestService.KEGGgetMol(candidate);
			MDLReader reader;
			List<IAtomContainer> containersList;
			
	        reader = new MDLReader(new StringReader(candidateMol));
	        ChemFile chemFile = (ChemFile)reader.read((ChemObject)new ChemFile());
	        containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
	        molecule = containersList.get(0);
			
		}
		else if(database.equals("chemspider"))
		{
			String candidateMol = ChemSpider.getMolByID(candidate, chemspiderToken);
			
			MDLReader reader;
			List<IAtomContainer> containersList;
			
	        reader = new MDLReader(new StringReader(candidateMol));
	        ChemFile chemFile = (ChemFile)reader.read((ChemObject)new ChemFile());
	        containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
	        molecule = containersList.get(0);
	        
		}
		else if(database.equals("pubchem"))
		{
			molecule = pw.getMol(candidate);
		}
		else
		{
			System.err.println("No database selected or wrong database name?");
		}
		
		return molecule;
	}
	
	/**
	 * Gets the compound locally.
	 * 
	 * @param database the database
	 * @param candidate the candidate
	 * @param jdbc the jdbc
	 * @param username the username
	 * @param password the password
	 * @param getAll the get all
	 * 
	 * @return the compound locally
	 * 
	 * @throws ClassNotFoundException the class not found exception
	 * @throws SQLException the SQL exception
	 * @throws CDKException the CDK exception
	 * @throws RemoteException the remote exception
	 */
	public static IAtomContainer getCompoundLocally(String database, String candidate, String jdbc, String username, String password, 
			boolean getAll, String chemspiderToken) throws SQLException, ClassNotFoundException, RemoteException, CDKException
	{
		IAtomContainer molecule = null;

		if(database.equals("kegg"))
		{
			KEGGLocal kl = new KEGGLocal(jdbc, username, password);
			molecule = kl.getMol(candidate);
		}
		else if(database.equals("chemspider"))
		{
			molecule = ChemSpider.getMol(candidate, chemspiderToken, getAll);
		}
		else if(database.equals("pubchem"))
		{
			PubChemLocal pl = new PubChemLocal(jdbc, username, password);
			molecule = pl.getMol(candidate, getAll);
		}
		
		return molecule;
	}
	
	/**
	 * removes structures whose first inchi key strings are equal
	 * 
	 * @return
	 * @throws CDKException 
	 * @throws RemoteException 
	 */
	public static Vector<String> removeDuplicatesByInchi(Vector<String> candidates, String database, 
			String chemspiderToken) throws RemoteException, CDKException {
		
		boolean[] uniqueStructures = new boolean[candidates.size()];
		for(int i = 0; i < uniqueStructures.length; i++) uniqueStructures[i] = true;
		
		//for inchi calculation the IAtomContainers are needed
		IAtomContainer[] molecules = new IAtomContainer[candidates.size()];
		
		for(int i = 0; i < candidates.size(); i++)
			molecules[i] = Candidates.getCompound(database, candidates.get(i), null, chemspiderToken);

		InChIGeneratorFactory igf = InChIGeneratorFactory.getInstance();
		for(int index1 = 0; index1 < molecules.length; index1++) {
			
			if(uniqueStructures[index1]) {
				
				InChIGenerator ig = igf.getInChIGenerator(molecules[index1]);
				if(ig.getReturnStatus() == INCHI_RET.ERROR) {
					System.err.println("Error creating InChI for [" + candidates.get(index1) + "]");
					continue;
				}
				
				String inchikey1 = ig.getInchiKey();
				if(inchikey1 == null || inchikey1.isEmpty())
					continue;
				String inchi1 = inchikey1.split("-")[0];
				
				for(int index2 = index1 + 1; index2 < molecules.length; index2++) {
					
					if(uniqueStructures[index2]) {
						
						ig = igf.getInChIGenerator(molecules[index2]);
						String inchikey2 = ig.getInchiKey();
						if(inchikey2 == null || inchikey2.isEmpty())
							continue;
						String inchi2 = inchikey2.split("-")[0];
					
						if(inchi1.compareTo(inchi2) == 0) {
							uniqueStructures[index2] = false;
						}
					
					
					}
				}
			}
		}
		
		Vector<String> cleanedCandidates = new Vector<String>();
		for(int i = 0; i < uniqueStructures.length; i++) 
			if(uniqueStructures[i]) cleanedCandidates.add(candidates.get(i));
		
		return cleanedCandidates;
	
	}
	
	/**
	 * removes structures whose first inchi key strings are equal
	 * 
	 * @return
	 * @throws CDKException 
	 * @throws RemoteException 
	 * @throws ClassNotFoundException 
	 * @throws SQLException 
	 */
	public static Vector<String> removeDuplicatesByInchiLocally(List<String> candidates, String database, String databaseUrl, 
			String username, String password, String chemspiderToken) throws RemoteException, CDKException, SQLException, ClassNotFoundException {
		
		boolean[] uniqueStructures = new boolean[candidates.size()];
		for(int i = 0; i < uniqueStructures.length; i++) uniqueStructures[i] = true;
		
		//for inchi calculation the IAtomContainers are needed
		IAtomContainer[] molecules = new IAtomContainer[candidates.size()];
		
		for(int i = 0; i < candidates.size(); i++)
			molecules[i] = Candidates.getCompoundLocally(database, candidates.get(i), databaseUrl, username, password, false, chemspiderToken);

		InChIGeneratorFactory igf = InChIGeneratorFactory.getInstance();
		for(int index1 = 0; index1 < molecules.length; index1++) {
			
			if(uniqueStructures[index1]) {
				
				InChIGenerator ig = igf.getInChIGenerator(molecules[index1]);
				if(ig.getReturnStatus() == INCHI_RET.ERROR) {
					System.err.println("Error creating InChI for [" + candidates.get(index1) + "]");
					continue;
				}

				String inchikey1 = ig.getInchiKey();
				if(inchikey1 == null || inchikey1.isEmpty())
					continue;
				
				String inchi1 = inchikey1.split("-")[0];
				
				for(int index2 = index1 + 1; index2 < molecules.length; index2++) {
					
					if(uniqueStructures[index2]) {
						
						ig = igf.getInChIGenerator(molecules[index2]);
						String inchikey2 = ig.getInchiKey();
						if(inchikey2 == null || inchikey2.isEmpty())
							continue;
						String inchi2 = inchikey2.split("-")[0];
					
						if(inchi1.compareTo(inchi2) == 0) {
							uniqueStructures[index2] = false;
						}
					
					
					}
				}
			}
		}
		
		Vector<String> cleanedCandidates = new Vector<String>();
		for(int i = 0; i < uniqueStructures.length; i++) 
			if(uniqueStructures[i]) cleanedCandidates.add(candidates.get(i));
		
		return cleanedCandidates;
	
	}
}
