
/**
 * PUGSkeleton.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis2 version: 1.4.1  Built on : Aug 13, 2008 (05:03:35 LKT)
 */
    package gov.nih.nlm.ncbi.pubchemAxis2;
    /**
     *  PUGSkeleton java skeleton for the axisService
     */
    public class PUGSkeleton{
        
         
        /**
         * Auto generated method signature
         * Given a list key, return the identifiers as an array of integers. Note that this method expects there to be at least one identifier in the list, and will fault if the list is empty; see GetListItemsCount, which can be used to check for an empty list prior to calling GetIDList. The optional Start (zero-based) and Count parameters can be used to return smaller portions of the list, useful especially for large lists. Synchronous.
                                     * @param getIDList
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetIDListResponse GetIDList
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetIDList getIDList
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetIDList");
        }
     
         
        /**
         * Auto generated method signature
         * Get the description of column (readout) in a BioAssay, which may be the outcome, score, or a TID from the given AID. Synchronous.
                                     * @param getAssayColumnDescription
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetAssayColumnDescriptionResponse GetAssayColumnDescription
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetAssayColumnDescription getAssayColumnDescription
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetAssayColumnDescription");
        }
     
         
        /**
         * Auto generated method signature
         * Get the descriptive information for a BioAssay, including the number of user-specified readouts (TIDs) and whether a score readout is present. Optionally get version and SID/CID count information. If GetFullDataBlob is set to true, then a Base64 string of data is returned in the response instead, containing the full PubChem Assay description in the requested format (ASN or XML only). Synchronous.
                                     * @param getAssayDescription
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetAssayDescriptionResponse GetAssayDescription
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetAssayDescription getAssayDescription
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetAssayDescription");
        }
     
         
        /**
         * Auto generated method signature
         * Input a chemical structure in ASN.1 (text or binary), XML, or SDF format. The structure must be encoded as a Base64 string. Currently only single structures are supported. Returns a structure key. Synchronous.
                                     * @param inputStructureBase64
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputStructureBase64Response InputStructureBase64
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputStructureBase64 inputStructureBase64
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputStructureBase64");
        }
     
         
        /**
         * Auto generated method signature
         * Search PubChem Compound for structures containing the one given by the structure key input, based on a user-selected level of chemical identity: connectivity only, match isotopes and/or stereo, etc. The search may be limited by elapsed time or number of records found, or restricted to search only within a previous result set (given by a list key). Returns a list key. Asynchronous.
                                     * @param substructureSearch
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.SubstructureSearchResponse SubstructureSearch
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.SubstructureSearch substructureSearch
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#SubstructureSearch");
        }
     
         
        /**
         * Auto generated method signature
         * Search PubChem Compound for structures similar to the one given by the structure key input, based on the given Tanimoto-based similarity score. The search may be limited by elapsed time or number of records found, or restricted to search only within a previous result set (given by a list key). Returns a list key. Asynchronous.
                                     * @param similaritySearch2D
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.SimilaritySearch2DResponse SimilaritySearch2D
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.SimilaritySearch2D similaritySearch2D
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#SimilaritySearch2D");
        }
     
         
        /**
         * Auto generated method signature
         * Given a structure key that has been processed by Standardize, return the chemical structure as ASN, XML, or SDF, returned as a Base64-encoded string. Synchronous.
                                     * @param getStandardizedStructureBase64
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedStructureBase64Response GetStandardizedStructureBase64
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedStructureBase64 getStandardizedStructureBase64
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetStandardizedStructureBase64");
        }
     
         
        /**
         * Auto generated method signature
         * Search PubChem Compound for structures identical to the one given by the structure key input, based on a user-selected level of chemical identity: connectivity only, match isotopes and/or stereo, etc. The search may be limited by elapsed time or number of records found, or restricted to search only within a previous result set (given by a list key). Returns a list key. Asynchronous.
                                     * @param identitySearch
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.IdentitySearchResponse IdentitySearch
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.IdentitySearch identitySearch
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#IdentitySearch");
        }
     
         
        /**
         * Auto generated method signature
         * Return the number of IDs in the set represented by the given list key. Synchronous.
                                     * @param getListItemsCount
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetListItemsCountResponse GetListItemsCount
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetListItemsCount getListItemsCount
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetListItemsCount");
        }
     
         
        /**
         * Auto generated method signature
         * Given a list key, return an Entrez history key (db, query key, and WebEnv) corresponding to that list. Synchronous.
                                     * @param getEntrezKey
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetEntrezKeyResponse GetEntrezKey
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetEntrezKey getEntrezKey
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetEntrezKey");
        }
     
         
        /**
         * Auto generated method signature
         * Specify an assay table from a BioAssay AID. The table may be complete, concise, or include a ListKey-specified set of readouts (TIDs). By default, all tested substances are included, but can be restricted to a ListKey-specified set of SIDs or CIDs. Returns an assay key. Synchronous.
                                     * @param inputAssay
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputAssayResponse InputAssay
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputAssay inputAssay
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputAssay");
        }
     
         
        /**
         * Auto generated method signature
         * Input a set of identifiers for a PubChem database, as an array of integers. Returns a list key. Synchronous.
                                     * @param inputList
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputListResponse InputList
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputList inputList
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputList");
        }
     
         
        /**
         * Auto generated method signature
         * Search PubChem Compound for structures of a given molecular formula, optionally allowing elements not specified to be present. The search may be limited by elapsed time or number of records found, or restricted to search only within a previous result set (given by a list key). Returns a list key. Asynchronous.
                                     * @param mFSearch
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.MFSearchResponse MFSearch
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.MFSearch mFSearch
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#MFSearch");
        }
     
         
        /**
         * Auto generated method signature
         * Given a structure key that has been processed by Standardize, return the chemical structure in as SMILES or InChI strings. Synchronous.
                                     * @param getStandardizedStructure
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedStructureResponse GetStandardizedStructure
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedStructure getStandardizedStructure
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetStandardizedStructure");
        }
     
         
        /**
         * Auto generated method signature
         * Given a key for any asynchronous operation, return any system messages (error messages, job info, etc.) associated with the operation, if any. Synchronous.
                                     * @param getStatusMessage
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetStatusMessageResponse GetStatusMessage
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetStatusMessage getStatusMessage
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetStatusMessage");
        }
     
         
        /**
         * Auto generated method signature
         * Input an Entrez history key (db, query key, and WebEnv). Returns a list key. Synchronous.
                                     * @param inputEntrez
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputEntrezResponse InputEntrez
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputEntrez inputEntrez
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputEntrez");
        }
     
         
        /**
         * Auto generated method signature
         * Given a list key, prepare for download a file containing those records in the selected format. See the web download service documentation (http://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch-help.html) for more detail on the supported formats and file types. Returns a download key. Asynchronous. Note that if SynchronousSingleRecord is set to true, and the ListKey contains only a single ID, then a Base64 string of data is returned synchronously in the response, instead of going through the download file.
                                     * @param download
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.DownloadResponse Download
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.Download download
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#Download");
        }
     
         
        /**
         * Auto generated method signature
         * Given a download key, return an FTP URL that may be used to download the requested file. Synchronous.
                                     * @param getDownloadUrl
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetDownloadUrlResponse GetDownloadUrl
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetDownloadUrl getDownloadUrl
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetDownloadUrl");
        }
     
         
        /**
         * Auto generated method signature
         * Input a chemical structure as a simple (one-line) string, either SMILES or InChI. Returns a structure key. Synchronous.
                                     * @param inputStructure
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputStructureResponse InputStructure
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputStructure inputStructure
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputStructure");
        }
     
         
        /**
         * Auto generated method signature
         * Given a key for any asynchronous operation, return the status of that operation. Possible return values are: Success, the operation completed normally; HitLimit, TimeLimit: the operation finished normally, but one of the limits was reached (e.g. before the entire database was searched); ServerError, InputError, DataError, Stopped: there was a problem with the input or on the server, and the job has died; Queued: the operation is waiting its turn in the public queue; Running: the operation is in progress. Synchronous.
                                     * @param getOperationStatus
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetOperationStatusResponse GetOperationStatus
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetOperationStatus getOperationStatus
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetOperationStatus");
        }
     
         
        /**
         * Auto generated method signature
         * Input a set of identifiers for a PubChem database, as a simple string of integer values separated by commas and/or whitespace. Returns a list key. Synchronous.
                                     * @param inputListText
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputListTextResponse InputListText
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputListText inputListText
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputListText");
        }
     
         
        /**
         * Auto generated method signature
         * Given an Entrez history key (db, query key, and WebEnv), return an HTTP URL that may be used to view the list in Entrez. Synchronous.
                                     * @param getEntrezUrl
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetEntrezUrlResponse GetEntrezUrl
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetEntrezUrl getEntrezUrl
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetEntrezUrl");
        }
     
         
        /**
         * Auto generated method signature
         * Input a set of identifiers for a PubChem database, as an array of strings. Returns a list key. Synchronous.
                                     * @param inputListString
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.InputListStringResponse InputListString
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.InputListString inputListString
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#InputListString");
        }
     
         
        /**
         * Auto generated method signature
         * Convert IDs from one type to another, using any one of a variety of CID matching algorithms. Output can be a list or a downloaded file; download file compression is optional and defaults to gzip (.gz). Returns a list or download key. Asynchronous.
                                     * @param iDExchange
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.IDExchangeResponse IDExchange
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.IDExchange iDExchange
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#IDExchange");
        }
     
         
        /**
         * Auto generated method signature
         * Standardize the structure given by the structure key input, using the same algorithm PubChem uses to construct the Compound database. Returns a structure key. Asynchronous.
                                     * @param standardize
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.StandardizeResponse Standardize
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.Standardize standardize
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#Standardize");
        }
     
         
        /**
         * Auto generated method signature
         * Get the description of all columns (readouts) in a BioAssay. Synchronous.
                                     * @param getAssayColumnDescriptions
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetAssayColumnDescriptionsResponse GetAssayColumnDescriptions
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetAssayColumnDescriptions getAssayColumnDescriptions
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetAssayColumnDescriptions");
        }
     
         
        /**
         * Auto generated method signature
         * Given a structure key that has been processed by Standardize, return the corresponding PubChem Compound database CID, or an empty value if the structure is not present in PubChem. Synchronous.
                                     * @param getStandardizedCID
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedCIDResponse GetStandardizedCID
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.GetStandardizedCID getStandardizedCID
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#GetStandardizedCID");
        }
     
         
        /**
         * Auto generated method signature
         * Given an assay key, prepare for download a file containing an assay data table in the selected format. See the assay query section of the PUG service documentation (http://pubchem.ncbi.nlm.nih.gov/pug/pughelp.html) for more detail on the supported formats. Compression is optional and defaults to gzip (.gz). Returns a download key. Asynchronous.
                                     * @param assayDownload
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.AssayDownloadResponse AssayDownload
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.AssayDownload assayDownload
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#AssayDownload");
        }
     
         
        /**
         * Auto generated method signature
         * Search PubChem Compound for structures contained within the one given by the structure key input, based on a user-selected level of chemical identity: connectivity only, match isotopes and/or stereo, etc. The search may be limited by elapsed time or number of records found, or restricted to search only within a previous result set (given by a list key). Returns a list key. Asynchronous.
                                     * @param superstructureSearch
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.SuperstructureSearchResponse SuperstructureSearch
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.SuperstructureSearch superstructureSearch
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#SuperstructureSearch");
        }
     
         
        /**
         * Auto generated method signature
         * Compute a matrix of scores from one or two lists of IDs (if one, the IDs will be self-scored), of the selected type and in the selected format. Compression is optional and defaults to gzip (.gz). Returns a download key. Asynchronous.
                                     * @param scoreMatrix
         */
        
                 public gov.nih.nlm.ncbi.pubchemAxis2.ScoreMatrixResponse ScoreMatrix
                  (
                  gov.nih.nlm.ncbi.pubchemAxis2.ScoreMatrix scoreMatrix
                  )
            {
                //TODO : fill this with the necessary business logic
                throw new  java.lang.UnsupportedOperationException("Please implement " + this.getClass().getName() + "#ScoreMatrix");
        }
     
    }
    