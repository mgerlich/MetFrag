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
package de.ipbhalle.metfrag.chemspiderClient;

import java.io.StringReader;
import java.rmi.RemoteException;
import java.util.List;
import java.util.Vector;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import com.chemspider.www.ExtendedCompoundInfo;
import com.chemspider.www.MassSpecAPISoapProxy;

import de.ipbhalle.metfrag.tools.MolecularFormulaTools;

public class ChemSpider {
	
	private static String token = "eeca1d0f-4c03-4d81-aa96-328cdccf171a"; 
		//"a1004d0f-9d37-47e0-acdd-35e58e34f603";
	
	/**
	 * Gets Chemspider compound ID's by mass.
	 * 
	 * @param mass the mass
	 * @param error the error
	 * 
	 * @return the chemspider by mass
	 * 
	 * @throws RemoteException the remote exception
	 */
	public static Vector<String> getChemspiderByMass(Double mass, Double error) throws RemoteException
	{
		Vector<String> result = new Vector<String>();
		MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
		String[] resultIDs = chemSpiderProxy.searchByMass2(mass, error);
		for (int i = 0; i < resultIDs.length; i++) {
			result.add(resultIDs[i]);
		}
		return result;
	}
	
	/**
	 * Gets the chemspider by sum formula.
	 * 
	 * @param sumFormula the sum formula
	 * 
	 * @return the chemspider by sum formula
	 * 
	 * @throws RemoteException the remote exception
	 */
	public static Vector<String> getChemspiderBySumFormula(String sumFormula) throws RemoteException
	{
		Vector<String> result = new Vector<String>();
		MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
		String[] resultIDs = chemSpiderProxy.searchByFormula2(sumFormula);
		for (int i = 0; i < resultIDs.length; i++) {
			result.add(resultIDs[i]);
		}
		return result;
	}
	
	
	/**
	 * Gets the mol by id.
	 * 
	 * @param ID the iD
	 * 
	 * @return the molby id
	 * 
	 * @throws RemoteException the remote exception
	 */
	public static String getMolByID(String ID) throws RemoteException
	{
		MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
		String mol = chemSpiderProxy.getRecordMol(ID, false, token);
		return mol;
	}
	
	/**
	 * Gets the mol by id.
	 * 
	 * @param ID the iD
	 * 
	 * @return the molby id
	 * 
	 * @throws RemoteException the remote exception
	 * @throws CDKException 
	 */
	public static IAtomContainer getMol(String ID, boolean getAll) throws RemoteException, CDKException
	{
		MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
		String mol = chemSpiderProxy.getRecordMol(ID, false, token);
		MDLReader reader;
		List<IAtomContainer> containersList;
		
        reader = new MDLReader(new StringReader(mol));
        ChemFile chemFile = (ChemFile)reader.read((ChemObject)new ChemFile());
        containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
        IAtomContainer molecule = containersList.get(0);
		
        if(getAll)
        	return molecule;
        
        if(!MolecularFormulaTools.isBiologicalCompound(molecule))
    		molecule = null;
        
		return molecule;
	}
	
	/**
	 * Gets the extended compound information like name, mass, InchI.....
	 * 
	 * @param key the key
	 * 
	 * @return the extended cpd info
	 * 
	 * @throws RemoteException the remote exception
	 */
	public static ExtendedCompoundInfo getExtendedCpdInfo(int key) throws RemoteException
	{
		MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
		ExtendedCompoundInfo cpdInfo = chemSpiderProxy.getExtendedCompoundInfo(key, token);
		return cpdInfo;
	}
	
	
	public static void main(String[] args) {
		try {
			MassSpecAPISoapProxy chemSpiderProxy = new MassSpecAPISoapProxy();
			String[] resultIDs = chemSpiderProxy.searchByMass2(272.06847, 0.01);
			
			ExtendedCompoundInfo cpdInfo = chemSpiderProxy.getExtendedCompoundInfo(907, "4d6c67db-65d0-474e-9f5c-f70f5c85111c");
			
			for (int i = 0; i < resultIDs.length; i++) {
				System.out.println("Result: " + resultIDs[i]);
			}
			System.out.println("#Hits: " + resultIDs.length);
			System.out.println("Name: " + cpdInfo.getCommonName());
			System.out.println("MW: " + cpdInfo.getMolecularWeight());
			System.out.println("MM: " + cpdInfo.getMonoisotopicMass());
			
			
		}
		catch (Exception e) {
			System.err.println("ERROR: " + e.getMessage());
		}			
	}
}

