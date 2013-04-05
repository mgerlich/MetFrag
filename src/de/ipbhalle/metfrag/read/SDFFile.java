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
package de.ipbhalle.metfrag.read;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import de.ipbhalle.metfrag.main.MetFrag;

public class SDFFile {
	
	/**
	 * Read sdf file.
	 * 
	 * @param path the path
	 * 
	 * @return the list< i atom container>
	 * 
	 * @throws FileNotFoundException the file not found exception
	 * @throws CDKException the CDK exception
	 */
	public static List<IAtomContainer> ReadSDFFile(String path) throws FileNotFoundException
	{
		//MDLV2000Reader reader;
		MDLReader reader;
		List<IAtomContainer> containersList;
		List<IAtomContainer> ret = new ArrayList<IAtomContainer>();
		
		File f = new File(path);
		
		if(f.isFile())
		{
			//reader = new MDLV2000Reader(new FileReader(f));
			reader = new MDLReader(new FileInputStream(f));
			ChemFile chemFile = null;
			try {
				chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			} catch (CDKException e) {
				System.err.println("Error reading SDF file " + f.getAbsolutePath());
				return ret;
			}
			containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
			for (IAtomContainer container : containersList) {
				ret.add(container);
			}
		}
		else
		{
			System.err.println("Did not find SDF file: " + path);
			//throw error
//			reader = new MDLV2000Reader(new FileReader(f));
		}
		
        return ret;
	}
	
	public static void main(String[] args) {
		try {
			List<IAtomContainer> list = ReadSDFFile("/vol/data_extern/emma.schymanski@ufz.de/ufzleipzig/100spec/001_61627_C9H16_struct_wM_END.sdf");
			for (IAtomContainer molecule : list) {
				try
		        {
			        //add hydrogens
			        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
			        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
			        hAdder.addImplicitHydrogens(molecule);
			        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
		        }
		        //there is a bug in cdk??
		        catch(IllegalArgumentException e)
	            {
		        	System.err.println(e.getMessage());
		        	e.printStackTrace();
	            }
			}
			
			System.out.println(list.size());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
