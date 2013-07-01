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

package de.ipbhalle.metfrag.pubchem;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import javax.xml.rpc.ServiceException;

import net.sf.jniinchi.INCHI_RET;

import org.apache.axis2.AxisFault;
import org.apache.axis2.transport.http.HTTPConstants;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import gov.nih.nlm.ncbi.pubchem.CompressType;
import gov.nih.nlm.ncbi.pubchem.EntrezKey;
import gov.nih.nlm.ncbi.pubchem.FormatType;
import gov.nih.nlm.ncbi.pubchem.MFSearchOptions;
import gov.nih.nlm.ncbi.pubchem.PCIDType;
import gov.nih.nlm.ncbi.pubchem.PUGLocator;
import gov.nih.nlm.ncbi.pubchem.PUGSoap;
import gov.nih.nlm.ncbi.pubchem.StatusType;
import gov.nih.nlm.ncbi.pubchemAxis2.AnyKeyType;
import gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfInt;
import gov.nih.nlm.ncbi.pubchemAxis2.Download;
import gov.nih.nlm.ncbi.pubchemAxis2.DownloadResponse;
import gov.nih.nlm.ncbi.pubchemAxis2.GetDownloadUrl;
import gov.nih.nlm.ncbi.pubchemAxis2.GetOperationStatus;
import gov.nih.nlm.ncbi.pubchemAxis2.InputList;
import gov.nih.nlm.ncbi.pubchemAxis2.InputListResponse;
import gov.nih.nlm.ncbi.pubchemAxis2.PUGStub;
import gov.nih.nlm.ncbi.www.soap.eutils.EUtilsServiceLocator;
import gov.nih.nlm.ncbi.www.soap.eutils.EUtilsServiceSoap;
import gov.nih.nlm.ncbi.www.soap.eutils.esearch.ESearchRequest;
import gov.nih.nlm.ncbi.www.soap.eutils.esearch.ESearchResult;

public class PubChemWebService {
	
	EUtilsServiceLocator eutils_locator;
	EUtilsServiceSoap eutils_soap;
	PUGLocator pug_locator;
	PUGSoap pug_soap;
	List<IAtomContainer> containers;
	HashMap<Integer, String> retrievedHits = null;
	boolean verbose = false;
	
	/**
	 * Instantiates a new pub chem web service.
	 * 
	 * @throws ServiceException the service exception
	 */
	public PubChemWebService() throws ServiceException
	{
		eutils_locator = new EUtilsServiceLocator();
		eutils_soap = eutils_locator.geteUtilsServiceSoap();
		pug_locator = new PUGLocator();
		pug_soap = pug_locator.getPUGSoap();
		this.retrievedHits = new HashMap<Integer, String>();
		this.containers = new ArrayList<IAtomContainer>();
	}
	
	
	
	public IAtomContainer getSingleMol(String cid, boolean useProxy) throws CDKException, InterruptedException, IOException
	{
		IAtomContainer ac = null;
		int[] cids = new int[1];
		cids[0] = Integer.parseInt(cid);

        String listKey = pug_soap.inputList(cids, PCIDType.eID_CID);
//        System.out.println("ListKey = " + listKey);
//        System.out.println("number of compounds = " + pug_soap.getListItemsCount(listKey));
        
        // Initialize the download; request SDF with gzip compression
        String downloadKey = pug_soap.download(listKey, 
            FormatType.eFormat_SDF, CompressType.eCompress_GZip, false);
//        System.out.println("DownloadKey = " + downloadKey);
        
        // Wait for the download to be prepared
        StatusType status;
        while ((status = pug_soap.getOperationStatus(downloadKey)) == StatusType.eStatus_Running || 
               status == StatusType.eStatus_Queued) 
        {
//            System.out.println("Waiting for download to finish...");
            Thread.sleep(5000);
        }
        
        // On success, get the download URL, save to local file
        if (status == StatusType.eStatus_Success) {
        	
        	// PROXY
        	if(useProxy)
        	{
        		System.getProperties().put( "ftp.proxySet", "true" );
    		    System.getProperties().put( "ftp.proxyHost", "www.ipb-halle.de" );
    		    System.getProperties().put( "ftp.proxyPort", "3128" );
        	}
		    
            URL url = new URL(pug_soap.getDownloadUrl(downloadKey));
//            System.out.println("Success! Download URL = " + url.toString());
            
            // get input stream from URL
            URLConnection fetch = url.openConnection();
            InputStream input = fetch.getInputStream();
            
            // open local file based on the URL file name
            File tempFile = File.createTempFile(url.getFile().substring(url.getFile().lastIndexOf(System.getProperty("file.separator"))), ".sdf");

            // Delete temp file when program exits.
            tempFile.deleteOnExit();
            FileOutputStream output = new FileOutputStream(tempFile);

			// buffered read/write
			byte[] buffer = new byte[10000];
			int n;
			while ((n = input.read(buffer)) > 0)
				output.write(buffer, 0, n);
			output.close();
			//now read in the file
			FileInputStream in = null;
	        in = new FileInputStream(tempFile);
	        GZIPInputStream gin = new GZIPInputStream(in);
	        
	        //IChemObjectReader cor = null;
	        //cor = new ReaderFactory().createReader(in);
	       
	        MDLV2000Reader reader = new MDLV2000Reader(gin);
	        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
	        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
	        ac = fileContents.getChemSequence(0).getChemModel(0).getMoleculeSet().getAtomContainer(0);
	        
	        // close streams
	        reader.close();
	        gin.close();
	        in.close();
        } else {
            System.out.println("Error: ");            
        }
        
        return ac;
	}
	
	
	/**
	 * PubChem get hits by sum formula. TODO: fix inefficient smiles generation
	 * 
	 * @param sumFormula the sum formula
	 * 
	 * @return the vector< string>
	 * 
	 * @throws Exception the exception
	 */
	public Vector<String> getHitsbySumFormula(String sumFormula, boolean useProxy, boolean uniqueInchi) throws Exception
	{
		Vector<String> candidatesString = new Vector<String>();        
		
		PUGLocator pug_locator = new PUGLocator();
		PUGSoap pug_soap = pug_locator.getPUGSoap();
		
		
		MFSearchOptions mf_options = new MFSearchOptions();
		mf_options.setAllowOtherElements(false);
		String listKey = pug_soap.MFSearch(sumFormula, mf_options, null);
		
		System.out.println("MFSearch " + sumFormula + " " + listKey);
		
		StatusType status;
		while ((status = pug_soap.getOperationStatus(listKey)) == StatusType.eStatus_Running
				|| status == StatusType.eStatus_Queued) {
			System.out.println("Waiting for query to finish...");
			Thread.sleep(10000);
		}
		
		int[] cids = null;
		//get cids
		try
		{
			cids = pug_soap.getIDList(listKey);
		}
		catch(RemoteException e)
		{
			System.err.println("Error: No hit!?" + e.getMessage());
			return candidatesString;
		}
		
		String listkey = pug_soap.inputList(cids, PCIDType.eID_CID);
		String downloadKey = pug_soap.download(listkey, FormatType.eFormat_SDF,
				CompressType.eCompress_None, false);
		System.out.println("DownloadKey = " + downloadKey);
		status = null;
		while ((status = pug_soap.getOperationStatus(downloadKey)) == StatusType.eStatus_Running
				|| status == StatusType.eStatus_Queued) {
			System.out.println("Waiting for download to finish...");
			Thread.sleep(1000);
		}
		
		// On success, get the download URL, save to local file
		if (status == StatusType.eStatus_Success) {
			
			// PROXY
			if(useProxy)
			{
				System.getProperties().put( "ftp.proxySet", "true" );
			    System.getProperties().put( "ftp.proxyHost", "www.ipb-halle.de" );
			    System.getProperties().put( "ftp.proxyPort", "3128" );
			}
		    
			URL url = new URL(pug_soap.getDownloadUrl(downloadKey));
			System.out.println("Success! Download URL = " + url.toString());

			// get input stream from URL
			URLConnection fetch = url.openConnection();
			InputStream input = fetch.getInputStream();

			// open local file based on the URL file name
            File tempFile = File.createTempFile(url.getFile().substring(url.getFile().lastIndexOf(System.getProperty("file.separator"))), ".sdf");

            // Delete temp file when program exits.
            tempFile.deleteOnExit();
            FileOutputStream output = new FileOutputStream(tempFile);

			System.out.println("Writing data to " + tempFile.getName());

			// buffered read/write
			byte[] buffer = new byte[10000];
			int n;
			while ((n = input.read(buffer)) > 0)
				output.write(buffer, 0, n);
			output.close();
			
			if(uniqueInchi) {
				candidatesString = filterDuplicatesByInChIKey(tempFile);
			}
			else {
				//read the file
				FileInputStream in = null;
		        in = new FileInputStream(tempFile);
		        
		        //IChemObjectReader cor = null;
		        //cor = new ReaderFactory().createReader(in);
		       
		        MDLV2000Reader reader = new MDLV2000Reader(in);
		        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
		        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");

		        // close streams
		        reader.close();
		        in.close();
		        
		        //ReaderFactory factory = new ReaderFactory();
		        //ISimpleChemObjectReader reader = factory.createReader(in);
		        //IChemFile content = (IChemFile)reader.read(new ChemFile());
		        
		        //IChemFile content = (IChemFile)cor.read(DefaultChemObjectBuilder.getInstance().newChemFile());
		        
		        System.out.println("Read the file");
		        this.containers = ChemFileManipulator.getAllAtomContainers(fileContents);
		        System.out.println("Got " + containers.size() + " atom containers");
		        
		        // Retrieve CIDs
		        SmilesGenerator generatorSmiles = new SmilesGenerator();
				for (int i = 0; i < cids.length; i++) {
					candidatesString.add(cids[i] + "");
					System.out.println(cids[i]);
					this.retrievedHits.put(cids[i], generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
				}
			}
		}
		 else {
			System.out.println("Error: " + pug_soap.getStatusMessage(downloadKey));
		}
		
		return candidatesString;
	}
	
	
	/**
	 * Search pubchem by exact mass.
	 * 
	 * @param mass the mass
	 * 
	 * @return the pubchem by mass
	 * 
	 * @throws Exception the exception
	 */
	public Vector<String> getHitsByMass(double mass, double error, Integer limit, boolean useProxy) throws Exception
	{
		Vector<String> pubchemCIDs = new Vector<String>();
		
		PUGSoap pug_soap = this.pug_locator.getPUGSoap();
		
        EUtilsServiceLocator eutils_locator = new EUtilsServiceLocator();
        EUtilsServiceSoap eutils_soap = eutils_locator.geteUtilsServiceSoap();


		
		//search "aspirin" in PubChem Compound
		ESearchRequest request = new ESearchRequest();
		String db = new String("pccompound");
		request.setDb(db);
		
		double min = mass - error;
		double max = mass + error;
		
		System.out.println("Min: " + min + " Max: " + max);
		request.setTerm(min + ":" + max + "[EMAS]");
		// create a history item, and don't return any actual ids in the
		// SOAP response
		request.setUsehistory("y");
		request.setRetMax(limit.toString());

		ESearchResult result = eutils_soap.run_eSearch(request);
		
		//String[] idList = result.getIdList();
		if (result.getQueryKey() == null || result.getQueryKey().length() == 0
				|| result.getWebEnv() == null
				|| result.getWebEnv().length() == 0) {
			throw new Exception("ESearch failed to return query_key and WebEnv");
		}
		System.out.println("ESearch returned " + result.getCount() + " hits");
		
		
		// give this Entrez History info to PUG SOAP
		EntrezKey entrezKey = new EntrezKey(db, result.getQueryKey(), result.getWebEnv());	
		String listKey = pug_soap.inputEntrez(entrezKey);
		
		System.out.println("ListKey = " + listKey);
		
		//int[] ids = pug_soap.getIDList(entrezKey.getKey());
		
		// Initialize the download; request SDF with gzip compression
		String downloadKey = pug_soap.download(listKey, FormatType.eFormat_SDF,
				CompressType.eCompress_None, false);
		System.out.println("DownloadKey = " + downloadKey);

		// Wait for the download to be prepared
		StatusType status;
		while ((status = pug_soap.getOperationStatus(downloadKey)) == StatusType.eStatus_Running
				|| status == StatusType.eStatus_Queued) {
			System.out.println("Waiting for download to finish...");
			Thread.sleep(10000);
		}

		// On success, get the download URL, save to local file
		if (status == StatusType.eStatus_Success) {
			
			// PROXY
			if(useProxy)
			{
				System.getProperties().put( "ftp.proxySet", "true" );
			    System.getProperties().put( "ftp.proxyHost", "www.ipb-halle.de" );
			    System.getProperties().put( "ftp.proxyPort", "3128" );
			}
		    
			URL url = new URL(pug_soap.getDownloadUrl(downloadKey));
			System.out.println("Success! Download URL = " + url.toString());

			// get input stream from URL
			URLConnection fetch = url.openConnection();
			InputStream input = fetch.getInputStream();

			// open local file based on the URL file name
            File tempFile = File.createTempFile(url.getFile().substring(url.getFile().lastIndexOf(System.getProperty("file.separator"))), ".sdf");
            // Delete temp file when program exits.
            tempFile.deleteOnExit();
            FileOutputStream output = new FileOutputStream(tempFile);
			System.out.println("Writing data to " + tempFile.getAbsolutePath() + tempFile.getName());

			// buffered read/write
			byte[] buffer = new byte[10000];
			int n;
			while ((n = input.read(buffer)) > 0)
				output.write(buffer, 0, n);
			output.close();
			//now read in the file
			FileInputStream in = null;
	        in = new FileInputStream(tempFile);
	        
	        //IChemObjectReader cor = null;
	        //cor = new ReaderFactory().createReader(in);
	       
	        MDLV2000Reader reader = new MDLV2000Reader(in);
	        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
	        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
	        
	        // close streams
	        reader.close();
	        in.close();
	        
	        SmilesGenerator generatorSmiles = new SmilesGenerator();
	        for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
				this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
				Map<Object, Object> properties = fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
		        pubchemCIDs.add((String) properties.get("PUBCHEM_COMPOUND_CID"));
		        System.out.println((String) properties.get("PUBCHEM_COMPOUND_CID"));
		        this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
			}

	        System.out.println("Read the file");
			
		} else {
			System.out.println("Error: "
					+ pug_soap.getStatusMessage(downloadKey));
		}
		
		return pubchemCIDs;
		
	}
	
	/**
	 * Search pubchem by exact mass.
	 * 
	 * @param mass the mass
	 * 
	 * @return the pubchem by mass
	 * 
	 * @throws Exception the exception
	 */
	public Vector<String> getHitsByMass(double mass, double error, Integer limit, boolean useProxy, boolean uniqueInchi) throws Exception
	{
		Vector<String> pubchemCIDs = new Vector<String>();
		
		PUGSoap pug_soap = this.pug_locator.getPUGSoap();
		
        EUtilsServiceLocator eutils_locator = new EUtilsServiceLocator();
        EUtilsServiceSoap eutils_soap = eutils_locator.geteUtilsServiceSoap();

		//search "aspirin" in PubChem Compound
		ESearchRequest request = new ESearchRequest();
		String db = new String("pccompound");
		request.setDb(db);
		
		double min = mass - error;
		double max = mass + error;
		
		System.out.println("Min: " + min + " Max: " + max);
		request.setTerm(min + ":" + max + "[EMAS]");
		// create a history item, and don't return any actual ids in the
		// SOAP response
		request.setUsehistory("y");
		request.setRetMax(limit.toString());

		ESearchResult result = eutils_soap.run_eSearch(request);
		
		//String[] idList = result.getIdList();
		if (result.getQueryKey() == null || result.getQueryKey().length() == 0
				|| result.getWebEnv() == null
				|| result.getWebEnv().length() == 0) {
			throw new Exception("ESearch failed to return query_key and WebEnv");
		}
		System.out.println("ESearch returned " + result.getCount() + " hits");
		
		
		// give this Entrez History info to PUG SOAP
		EntrezKey entrezKey = new EntrezKey(db, result.getQueryKey(), result.getWebEnv());	
		String listKey = pug_soap.inputEntrez(entrezKey);
		
		System.out.println("ListKey = " + listKey);
		
		//int[] ids = pug_soap.getIDList(entrezKey.getKey());
		
		// Initialize the download; request SDF with gzip compression
		String downloadKey = pug_soap.download(listKey, FormatType.eFormat_SDF,
				CompressType.eCompress_None, false);
		System.out.println("DownloadKey = " + downloadKey);

		// Wait for the download to be prepared
		StatusType status;
		while ((status = pug_soap.getOperationStatus(downloadKey)) == StatusType.eStatus_Running
				|| status == StatusType.eStatus_Queued) {
			System.out.println("Waiting for download to finish...");
			Thread.sleep(10000);
		}

		// On success, get the download URL, save to local file
		if (status == StatusType.eStatus_Success) {
			
			// PROXY
			if(useProxy)
			{
				System.getProperties().put( "ftp.proxySet", "true" );
			    System.getProperties().put( "ftp.proxyHost", "www.ipb-halle.de" );
			    System.getProperties().put( "ftp.proxyPort", "3128" );
			}
		    
			URL url = new URL(pug_soap.getDownloadUrl(downloadKey));
			System.out.println("Success! Download URL = " + url.toString());
			
			// get input stream from URL
			URLConnection fetch = url.openConnection();
			InputStream input = fetch.getInputStream();
			
			// open local file based on the URL file name
            File tempFile = File.createTempFile(url.getFile().substring(url.getFile().lastIndexOf('/')), ".sdf");
            // Delete temp file when program exits.
            tempFile.deleteOnExit();
            FileOutputStream output = new FileOutputStream(tempFile);
			System.out.println("Writing data to " + tempFile.getAbsolutePath() + tempFile.getName());

			// buffered read/write
			byte[] buffer = new byte[10000];
			int n;
			while ((n = input.read(buffer)) > 0)
				output.write(buffer, 0, n);
			output.close();
			
			//now read in the file
			FileInputStream in = new FileInputStream(tempFile);
	        
	        //IChemObjectReader cor = null;
	        //cor = new ReaderFactory().createReader(in);
	       
	        MDLV2000Reader reader = new MDLV2000Reader(in);
	        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
	        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
	        
	        // close streams
	        reader.close();
	        in.close();
	        
	        SmilesGenerator generatorSmiles = new SmilesGenerator();
	        Vector<String> cids_inchi_keys = new Vector<String>();
	        for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
	        	Map<Object, Object> properties = fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
        	    String pubchemCID = (String) properties.get("PUBCHEM_COMPOUND_CID");
	        	boolean insert = true;
	        	this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
        		this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
        		if(uniqueInchi) {
	        		IAtomContainer molecule1 = getMol(pubchemCID);
					InChIGeneratorFactory igf = InChIGeneratorFactory.getInstance();
					if(igf != null) {
						InChIGenerator ig = igf.getInChIGenerator(molecule1);
						if(ig != null) {
							String inchi = ig.getInchiKey();
							if(inchi != null) {
								String first_part = inchi.split("-")[0];
								if(cids_inchi_keys.contains(first_part)) insert = false;
								else cids_inchi_keys.add(first_part);
							}
						}
					}
				}
	        	if(!insert) {
	        	    this.containers.remove(this.containers.size() - 1);
	        	    this.retrievedHits.remove(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()));
	        	}
	        	else {
	        		pubchemCIDs.add(pubchemCID);
	        	    System.out.println(pubchemCID);
	        	}
	        }
	        
	        System.out.println("Read the file");
			
		} else {
			System.out.println("Error: "
					+ pug_soap.getStatusMessage(downloadKey));
		}
		
		return pubchemCIDs;
		
	}
	
	/**
	 * when ids are given
	 * 
	 * @param ids
	 * @return
	 */
	public Vector<String> getHitsByIDs(Vector<String> ids) {
		File tempFile = null;
		try {
			tempFile = File.createTempFile(getRandomString(20), ".sdf");
			tempFile.deleteOnExit();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem CID download!");
			return new Vector<String>();
		}
		try {
			boolean success = savingRetrievedHits(tempFile, ids);
			if(!success) return new Vector<String>(); 
		} catch (AxisFault e) {
			System.err.println("Error: Could not open result stream when using Pubchem CID download!");
			return new Vector<String>();
		}
		return ids;
	}
	
	/**
	 * download directly over http because of server errors when downloading over pug soap
	 * 
	 * @author c-ruttkies
	 * 
	 * @param mass
	 * @param error
	 * @param limit
	 * @return
	 */
	public Vector<String> getHitsByMassHTTP(double mass, double error, int limit) {
		Vector<String> cids = new Vector<String>();
		double minMass = mass - error;
		double maxMass = mass + error;

		String urlname = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
				+ "?db=pccompound"
				+ "&term=" + minMass + "[ExactMass]:" + maxMass + "[ExactMass]"
				+ "&RetMax=" + limit;
		InputStream stream = getInputStreamFromURL(urlname);
		
		if(stream == null) return cids;
		try {
			BufferedReader breader = new BufferedReader(new InputStreamReader(stream));
			String line = "";
			while((line = breader.readLine()) != null) {
				if(line.contains("<Id>") && line.contains("</Id>")) {
					cids.add(line.replaceAll("\\D", "").trim());
				}
			}
			stream.close();
			breader.close();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
			return cids;
		}
		
		File tempFile = null;
		try {
			tempFile = File.createTempFile(getRandomString(20), ".sdf");
	        if(tempFile != null) tempFile.deleteOnExit();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
			return cids;
		}

		try {
			boolean success = savingRetrievedHits(tempFile, cids);
			if(!success) return new Vector<String>();
		} catch (AxisFault e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
		}
		
		return cids;
	}
	
	/**
	 * download directly over http because of server errors when downloading over pug soap
	 * 
	 * @author c-ruttkies
	 * 
	 * @param mass
	 * @param error
	 * @param limit
	 * @return
	 */
	public Vector<String> getHitsByMassHTTP(double mass, double error, int limit, boolean uniqueInchi) {
		Vector<String> cids = new Vector<String>();
		double minMass = mass - error;
		double maxMass = mass + error;

		String urlname = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
				+ "?db=pccompound"
				+ "&term=" + minMass + "[ExactMass]:" + maxMass + "[ExactMass]"
				+ "&RetMax=" + limit;
		InputStream stream = getInputStreamFromURL(urlname);
		
		if(stream == null) return cids;
		try {
			BufferedReader breader = new BufferedReader(new InputStreamReader(stream));
			String line = "";
			while((line = breader.readLine()) != null) {
				if(line.contains("<Id>") && line.contains("</Id>")) {
					cids.add(line.replaceAll("\\D", "").trim());
				}
			}
			stream.close();
			breader.close();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
			return cids;
		}
		
		File tempFile = null;
		try {
			tempFile = File.createTempFile(getRandomString(20), ".sdf");
	        if(tempFile != null) tempFile.deleteOnExit();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
			return cids;
		}

		try {
			boolean success = savingRetrievedHits(tempFile, cids);
			if(!success) return new Vector<String>();
		} catch (AxisFault e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
		}
		
		if(uniqueInchi) {
			cids = filterDuplicatesByInChIKey(tempFile);
			
//			cids = new Vector<String>();
//			
//			//now read in the file
//			FileInputStream in = null;
//	        try {
//				in = new FileInputStream(tempFile);
//				//IChemObjectReader cor = null;
//		        //cor = new ReaderFactory().createReader(in);
//		       
//		        MDLV2000Reader reader = new MDLV2000Reader(in);
//		        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
//		        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
//		        
//		        SmilesGenerator generatorSmiles = new SmilesGenerator();
//		        InChIGeneratorFactory igf = InChIGeneratorFactory.getInstance();
//		        Vector<String> cids_inchi_keys = new Vector<String>();
//		        for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
//		        	Map<Object, Object> properties = fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
//	        	    String pubchemCID = (String) properties.get("PUBCHEM_COMPOUND_CID");
//		        	boolean insert = true;
//		        	this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
//	        		this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
//	        		IAtomContainer molecule1 = getMol(pubchemCID);
//					if(igf != null) {
//						InChIGenerator ig = igf.getInChIGenerator(molecule1);
//						if(ig != null) {
//							if(ig.getReturnStatus().equals(INCHI_RET.OKAY) | ig.getReturnStatus().equals(INCHI_RET.WARNING)) {
//								String inchi = ig.getInchiKey();
//								if(inchi != null) {
//									String first_part = inchi.split("-")[0];
//									if(cids_inchi_keys.contains(first_part)) insert = false;
//									else cids_inchi_keys.add(first_part);
//								}
//							}
//						}
//					}
//					
//		        	if(!insert) {
//		        	    this.containers.remove(this.containers.size() - 1);
//		        	    this.retrievedHits.remove(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()));
//		        	    System.out.println("Removed " + pubchemCID);
//		        	}
//		        	else {
//		        		cids.add(pubchemCID);
//		        	    System.out.println(pubchemCID);
//		        	}
//		        }
//
//		        System.out.println("Read the file");
//			} catch (FileNotFoundException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			} catch (CDKException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
		}
		
		return cids;
	}
	
	private Vector<String> filterDuplicatesByInChIKey(File file) {
		Vector<String> cids = new Vector<String>();
		
		//now read in the file
		FileInputStream in = null;
        try {
			in = new FileInputStream(file);
			//IChemObjectReader cor = null;
	        //cor = new ReaderFactory().createReader(in);
	       
	        MDLV2000Reader reader = new MDLV2000Reader(in);
	        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
	        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
	        
	        // close streams
	        reader.close();
	        in.close();
	        
	        SmilesGenerator generatorSmiles = new SmilesGenerator();
	        InChIGeneratorFactory igf = InChIGeneratorFactory.getInstance();
	        Vector<String> cids_inchi_keys = new Vector<String>();
	        //for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
	        for (int i = fileContents.getChemSequence(0).getChemModelCount()-1; i >= 0 ; i--) {		// read in file in reverse order, smallest ID first
	        	Map<Object, Object> properties = fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
        	    String pubchemCID = (String) properties.get("PUBCHEM_COMPOUND_CID");
	        	boolean insert = true;
	        	this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
        		this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
        		IAtomContainer molecule1 = getMol(pubchemCID);
				if(igf != null) {
					InChIGenerator ig = igf.getInChIGenerator(molecule1);
					if(ig != null) {
						if(ig.getReturnStatus().equals(INCHI_RET.OKAY) | ig.getReturnStatus().equals(INCHI_RET.WARNING)) {
							String inchi = ig.getInchiKey();
							if(inchi != null) {
								String first_part = inchi.split("-")[0];
								if(cids_inchi_keys.contains(first_part)) insert = false;
								else cids_inchi_keys.add(first_part);
							}
						}
					}
				}
				
	        	if(!insert) {
	        	    this.containers.remove(this.containers.size() - 1);
	        	    this.retrievedHits.remove(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()));
	        	    System.out.println("Removed " + pubchemCID);
	        	}
	        	else {
	        		cids.add(pubchemCID);
	        	    System.out.println(pubchemCID);
	        	}
	        }

	        System.out.println("Read the file");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        return cids;
	}
	
	/**
	 * @author c-ruttkies
	 * 
	 * @param formula
	 * @return
	 */
	public Vector<String> getHitsBySumFormulaHTTP(String formula) {
		Vector<String> cids = new Vector<String>();
		String urlname = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
				+ "?db=pccompound"
				+ "&term=" + formula
				+ "&RetMax=10000000";
		InputStream stream = getInputStreamFromURL(urlname);
		if(stream == null) return cids;
		try {
			BufferedReader breader = new BufferedReader(new InputStreamReader(stream));
			String line = "";
			while((line = breader.readLine()) != null) {
				if(line.contains("<Id>") && line.contains("</Id>")) {
					cids.add(line.replaceAll("\\D", "").trim());
				}
			}
			stream.close();
			breader.close();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP formula search!");
			System.err.println(urlname);
			e.printStackTrace();
		}
		
		File tempFile = null;
		try {
			tempFile = File.createTempFile(getRandomString(20), ".sdf");
	        tempFile.deleteOnExit();
		} catch (IOException e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
		}

		try {
			boolean success = savingRetrievedHits(tempFile, cids);
			if(!success) return new Vector<String>();
		} catch (AxisFault e) {
			System.err.println("Error: Could not open result stream when using Pubchem HTTP mass search!");
			System.err.println(urlname);
		}

		return cids;
	}
	
	/**
	 * 
	 * @author c-ruttkies
	 * 
	 * 
	 * @param filename
	 * @throws AxisFault 
	 */
	private boolean savingRetrievedHits(File filename, Vector<String> cidsVec) throws AxisFault {
		org.apache.commons.httpclient.params.DefaultHttpParams.getDefaultParams().setParameter("http.protocol.cookie-policy",
				org.apache.commons.httpclient.cookie.CookiePolicy.BROWSER_COMPATIBILITY);

		
		PUGStub ps = new PUGStub();
		ps._getServiceClient().getOptions().setProperty(HTTPConstants.CHUNKED, false);
		ps._getServiceClient().getOptions().setTimeOutInMilliSeconds(5*60*1000);
		
		ArrayOfInt aoi = new ArrayOfInt();
		InputList il = new InputList();
		il.setIdType(gov.nih.nlm.ncbi.pubchemAxis2.PCIDType.eID_CID);
		Download d = new Download();
		d.setEFormat(gov.nih.nlm.ncbi.pubchemAxis2.FormatType.eFormat_SDF);
		d.setECompress(gov.nih.nlm.ncbi.pubchemAxis2.CompressType.eCompress_None);
		d.setUse3D(false);
		int[] cids = new int[cidsVec.size()];
		//for(int i = 0; i < cidsVec.size(); i++) {
		for(int i = cidsVec.size()-1; i >= 0 ; i--) {	// reverse order, lowest ID now first
			try {
				cids[i] = Integer.parseInt(cidsVec.get(i));
			} catch(java.lang.NumberFormatException e) {
				System.err.println("Error: "+cidsVec.get(i)+" is no valid pubchem ID!");
				return false;
			}
		}
       	aoi.set_int(cids);
		il.setIds(aoi);
		InputListResponse ilr = null;
		try {
			ilr = ps.InputList(il);
		} catch (RemoteException e) {
			System.err.println("Error: Pubchem sdf download failed. Contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false;
		} 
		if(ilr == null){
			System.err.println("Error: Pubchem sdf download failed. Contact cruttkies@ipb-halle.de!");
			return false;
		}
		String listKey = ilr.getListKey();
		d.setListKey(listKey);
		DownloadResponse dr = null;
		try {
			dr = ps.Download(d);
		} catch (RemoteException e) {
			System.err.println("Error: Pubchem sdf download failed. Contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false;
		}
		gov.nih.nlm.ncbi.pubchemAxis2.GetOperationStatus req4 = new GetOperationStatus();
		AnyKeyType anyKey = new AnyKeyType();
        anyKey.setAnyKey(dr.getDownloadKey());
        req4.setGetOperationStatus(anyKey);
        gov.nih.nlm.ncbi.pubchemAxis2.StatusType status;
        try {
			if(this.verbose) System.out.print("downloading compounds from pubchem");
			while ((status = ps.GetOperationStatus(req4).getStatus()) 
			        == gov.nih.nlm.ncbi.pubchemAxis2.StatusType.eStatus_Running || 
			   status == gov.nih.nlm.ncbi.pubchemAxis2.StatusType.eStatus_Queued) 
			{
				Thread.sleep(2000);
				if(this.verbose) System.out.print(".");
			}
		} catch (RemoteException e) {
			System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false;
		} catch (InterruptedException e) {
			System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false;
		}
        
        if (status == gov.nih.nlm.ncbi.pubchemAxis2.StatusType.eStatus_Success) {
            GetDownloadUrl req5 = new GetDownloadUrl();
            req5.setDownloadKey(dr.getDownloadKey());
            URL url = null;
			try {
				url = new URL(ps.GetDownloadUrl(req5).getUrl());
			} catch (MalformedURLException e) {
				System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
				e.printStackTrace();
				return false;
			} catch (RemoteException e) {
				System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
				e.printStackTrace();
				return false;
			}
            
			if(this.verbose) System.out.println("\ndownload finished!");
            
			URLConnection fetch;
            InputStream input;
			try {
				fetch = url.openConnection();
	            input = fetch.getInputStream();
			} catch (IOException e) {
				System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
				e.printStackTrace();
				return false;
			}
            
            FileOutputStream output;
			try {
				output = new FileOutputStream(filename.getAbsoluteFile());
			} catch (FileNotFoundException e) {
				System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
				e.printStackTrace();
				return false;
			}
            
            byte[] buffer = new byte[10000];
            int n;
            try {
				while ((n = input.read(buffer)) > 0)
				    output.write(buffer, 0, n);
				output.close();
			} catch (IOException e) {
				System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
				e.printStackTrace();
				return false;  
			}
        } else {
            System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
            return false;            
        }
        
        FileInputStream in = null;
        try {
			in = new FileInputStream(filename);
		} catch (FileNotFoundException e) {
			System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
			return false; 
		}
        
        MDLV2000Reader reader = new MDLV2000Reader(in);
        ChemFile fileContents = null;
		try {
			fileContents = (ChemFile)reader.read(new ChemFile());
			in.close();
			reader.close();
		} catch (CDKException e) {
			System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false; 
		} catch (IOException e) {
			System.err.println("Error: Pubchem sdf download failed. Please contact cruttkies@ipb-halle.de!");
			e.printStackTrace();
			return false; 
		}
        
        
        SmilesGenerator generatorSmiles = new SmilesGenerator();
        for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
			this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
			
			Map<Object, Object> properties = 
					fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
	        		this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), 
	        		generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
		}
        
        return true;
	}
	
	/**
	 * 
	 * @author c-ruttkies
	 * 
	 * @param urlname
	 * @return
	 */
	private static InputStream getInputStreamFromURL(String urlname) {
		InputStream stream = null;
		
		try {
			URL url = new URL(urlname);
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
		
			if (conn.getResponseCode() != 200) {
				throw new IOException(conn.getResponseMessage());
			}
			stream = conn.getInputStream();
		} catch(MalformedURLException mue) {
			System.err.println("Error: Could create URL object!");
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Error: Could not open URL connection!");
			System.exit(2);
		}
		
		return stream;
	}
	
	/**
	 * Gets pubchem candidate compounds by molecular weight.
	 * 
	 * @param mass the mass
	 * @param error the error
	 * @param limit the limit
	 * 
	 * @return the pubchem mimw
	 * 
	 * @throws Exception the exception
	 */
	public Vector<String> getHitsByMW(double mass, double error, Integer limit) throws Exception
	{
		Vector<String> pubchemCIDs = new Vector<String>();
		
		PUGSoap pug_soap = this.pug_locator.getPUGSoap();
		
        EUtilsServiceLocator eutils_locator = new EUtilsServiceLocator();
        EUtilsServiceSoap eutils_soap = eutils_locator.geteUtilsServiceSoap();


		
		//search "aspirin" in PubChem Compound
		ESearchRequest request = new ESearchRequest();
		String db = new String("pccompound");
		request.setDb(db);
		
		double min = mass - error;
		double max = mass + error;
		
		System.out.println("Min: " + min + " Max: " + max);
		
		request.setTerm(min + ":" + max + "[MW]");
		// create a history item, and don't return any actual ids in the
		// SOAP response
		request.setUsehistory("y");
		request.setRetMax(limit.toString());

		ESearchResult result = eutils_soap.run_eSearch(request);
		
		//String[] idList = result.getIdList();
		if (result.getQueryKey() == null || result.getQueryKey().length() == 0
				|| result.getWebEnv() == null
				|| result.getWebEnv().length() == 0) {
			throw new Exception("ESearch failed to return query_key and WebEnv");
		}
		System.out.println("ESearch returned " + result.getCount() + " hits");
		
		
		// give this Entrez History info to PUG SOAP
		EntrezKey entrezKey = new EntrezKey(db, result.getQueryKey(), result.getWebEnv());	
		String listKey = pug_soap.inputEntrez(entrezKey);
		
		System.out.println("ListKey = " + listKey);
		
		//int[] ids = pug_soap.getIDList(entrezKey.getKey());
		
		// Initialize the download; request SDF with gzip compression
		String downloadKey = pug_soap.download(listKey, FormatType.eFormat_SDF,
				CompressType.eCompress_None, false);
		System.out.println("DownloadKey = " + downloadKey);

		// Wait for the download to be prepared
		StatusType status;
		while ((status = pug_soap.getOperationStatus(downloadKey)) == StatusType.eStatus_Running
				|| status == StatusType.eStatus_Queued) {
			System.out.println("Waiting for download to finish...");
			Thread.sleep(10000);
		}

		// On success, get the download URL, save to local file
		if (status == StatusType.eStatus_Success) {
			
			// PROXY
		    System.getProperties().put( "ftp.proxySet", "true" );
		    System.getProperties().put( "ftp.proxyHost", "www.ipb-halle.de" );
		    System.getProperties().put( "ftp.proxyPort", "3128" );

			
			URL url = new URL(pug_soap.getDownloadUrl(downloadKey));
			System.out.println("Success! Download URL = " + url.toString());

			// get input stream from URL
			URLConnection fetch = url.openConnection();
			InputStream input = fetch.getInputStream();

			// open local file based on the URL file name
            File tempFile = File.createTempFile(url.getFile().substring(url.getFile().lastIndexOf('/')), ".sdf");
            // Delete temp file when program exits.
            tempFile.deleteOnExit();
            FileOutputStream output = new FileOutputStream(tempFile);
			System.out.println("Writing data to " + tempFile.getAbsolutePath() + tempFile.getName());

			// buffered read/write
			byte[] buffer = new byte[10000];
			int n;
			while ((n = input.read(buffer)) > 0)
				output.write(buffer, 0, n);
			output.close();
			//now read in the file
			FileInputStream in = null;
	        in = new FileInputStream(tempFile);
	        
	        //IChemObjectReader cor = null;
	        //cor = new ReaderFactory().createReader(in);
	       
	        MDLV2000Reader reader = new MDLV2000Reader(in);
	        ChemFile fileContents = (ChemFile)reader.read(new ChemFile());
	        System.out.println("Got " + fileContents.getChemSequence(0).getChemModelCount() + " atom containers");
	        
	        // close streams
	        reader.close();
	        in.close();
	        
	        SmilesGenerator generatorSmiles = new SmilesGenerator();
	        for (int i = 0; i < fileContents.getChemSequence(0).getChemModelCount(); i++) {
				this.containers.add(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0));
				Map<Object, Object> properties = fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getAtomContainer(0).getProperties();
		        pubchemCIDs.add((String) properties.get("PUBCHEM_COMPOUND_CID"));
		        System.out.println((String) properties.get("PUBCHEM_COMPOUND_CID"));
		        this.retrievedHits.put(Integer.parseInt(properties.get("PUBCHEM_COMPOUND_CID").toString()), generatorSmiles.createSMILES(fileContents.getChemSequence(0).getChemModel(i).getMoleculeSet().getMolecule(0)));
			}

	        System.out.println("Read the file");
			
		} else {
			System.out.println("Error: "
					+ pug_soap.getStatusMessage(downloadKey));
		}
		
		return pubchemCIDs;
		
	}
	
	
	
	/**
	 * Gets the compound. You have to execute the find by mass or find by
	 * molecular formula first, otherwise it will be null.
	 * 
	 * TODO
	 * 
	 * @param pubchemCID the pubchem cid
	 * 
	 * @return the compound
	 * @throws InvalidSmilesException 
	 */
	public IAtomContainer getCompound(int pubchemCID) throws InvalidSmilesException
	{
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		return sp.parseSmiles(this.retrievedHits.get(pubchemCID));		
	}
	
	
	/**
	 * Gets the mol. You have to execute the find by mass or find by
	 * molecular formula first, otherwise it will be null.
	 * 
	 * TODO
	 * 
	 * @param number the number
	 * 
	 * @return the mol
	 * @throws InvalidSmilesException 
	 * @throws NumberFormatException 
	 */
	public IAtomContainer getMol(String number) throws NumberFormatException, InvalidSmilesException
	{
		//got a new database hit...which is not stored in the database
		if(this.retrievedHits.get(Integer.parseInt(number)) == null)
			return null;
		
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		return sp.parseSmiles(this.retrievedHits.get(Integer.parseInt(number)));
	}

	/**
	 * 
	 * @param size
	 * @return
	 */
	private String getRandomString(int size) {
		char[] vals = {'0','1','2','3','4','5','6','7','8','9','Q','W','E','R','T','Z','U','I','O','P','A','S','D','F','G','H','J',
				'K','L','Y','X','C','V','B','N','M'};
		String randomString = "";
		java.util.Random rand = new java.util.Random();
		for(int i = 0; i < size; i++)
			randomString += vals[rand.nextInt(size)];
		
		return randomString;
	}
	
	
	public void setVerbose(boolean val) {
		this.verbose = val;
	}

}
