
/**
 * ColumnDescriptionType.java
 *
 * This file was auto-generated from WSDL
 * by the Apache Axis2 version: 1.4.1  Built on : Aug 13, 2008 (05:03:41 LKT)
 */
            
                package gov.nih.nlm.ncbi.pubchemAxis2;
            

            /**
            *  ColumnDescriptionType bean class
            */
        
        public  class ColumnDescriptionType
        implements org.apache.axis2.databinding.ADBBean{
        /* This type was generated from the piece of schema that had
                name = ColumnDescriptionType
                Namespace URI = http://pubchem.ncbi.nlm.nih.gov/
                Namespace Prefix = ns1
                */
            

        private static java.lang.String generatePrefix(java.lang.String namespace) {
            if(namespace.equals("http://pubchem.ncbi.nlm.nih.gov/")){
                return "ns1";
            }
            return org.apache.axis2.databinding.utils.BeanUtil.getUniquePrefix();
        }

        

                        /**
                        * field for Heading
                        */

                        
                                    protected gov.nih.nlm.ncbi.pubchemAxis2.HeadingType localHeading ;
                                

                           /**
                           * Auto generated getter method
                           * @return gov.nih.nlm.ncbi.pubchemAxis2.HeadingType
                           */
                           public  gov.nih.nlm.ncbi.pubchemAxis2.HeadingType getHeading(){
                               return localHeading;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param Heading
                               */
                               public void setHeading(gov.nih.nlm.ncbi.pubchemAxis2.HeadingType param){
                            
                                            this.localHeading=param;
                                    

                               }
                            

                        /**
                        * field for TID
                        */

                        
                                    protected int localTID ;
                                
                           /*  This tracker boolean wil be used to detect whether the user called the set method
                          *   for this attribute. It will be used to determine whether to include this field
                           *   in the serialized XML
                           */
                           protected boolean localTIDTracker = false ;
                           

                           /**
                           * Auto generated getter method
                           * @return int
                           */
                           public  int getTID(){
                               return localTID;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param TID
                               */
                               public void setTID(int param){
                            
                                       // setting primitive attribute tracker to true
                                       
                                               if (param==java.lang.Integer.MIN_VALUE) {
                                           localTIDTracker = false;
                                              
                                       } else {
                                          localTIDTracker = true;
                                       }
                                   
                                            this.localTID=param;
                                    

                               }
                            

                        /**
                        * field for Name
                        */

                        
                                    protected java.lang.String localName ;
                                

                           /**
                           * Auto generated getter method
                           * @return java.lang.String
                           */
                           public  java.lang.String getName(){
                               return localName;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param Name
                               */
                               public void setName(java.lang.String param){
                            
                                            this.localName=param;
                                    

                               }
                            

                        /**
                        * field for Description
                        */

                        
                                    protected gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfString localDescription ;
                                
                           /*  This tracker boolean wil be used to detect whether the user called the set method
                          *   for this attribute. It will be used to determine whether to include this field
                           *   in the serialized XML
                           */
                           protected boolean localDescriptionTracker = false ;
                           

                           /**
                           * Auto generated getter method
                           * @return gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfString
                           */
                           public  gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfString getDescription(){
                               return localDescription;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param Description
                               */
                               public void setDescription(gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfString param){
                            
                                       if (param != null){
                                          //update the setting tracker
                                          localDescriptionTracker = true;
                                       } else {
                                          localDescriptionTracker = false;
                                              
                                       }
                                   
                                            this.localDescription=param;
                                    

                               }
                            

                        /**
                        * field for Type
                        */

                        
                                    protected java.lang.String localType ;
                                

                           /**
                           * Auto generated getter method
                           * @return java.lang.String
                           */
                           public  java.lang.String getType(){
                               return localType;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param Type
                               */
                               public void setType(java.lang.String param){
                            
                                            this.localType=param;
                                    

                               }
                            

                        /**
                        * field for Unit
                        */

                        
                                    protected java.lang.String localUnit ;
                                
                           /*  This tracker boolean wil be used to detect whether the user called the set method
                          *   for this attribute. It will be used to determine whether to include this field
                           *   in the serialized XML
                           */
                           protected boolean localUnitTracker = false ;
                           

                           /**
                           * Auto generated getter method
                           * @return java.lang.String
                           */
                           public  java.lang.String getUnit(){
                               return localUnit;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param Unit
                               */
                               public void setUnit(java.lang.String param){
                            
                                       if (param != null){
                                          //update the setting tracker
                                          localUnitTracker = true;
                                       } else {
                                          localUnitTracker = false;
                                              
                                       }
                                   
                                            this.localUnit=param;
                                    

                               }
                            

                        /**
                        * field for TestedConcentration
                        */

                        
                                    protected gov.nih.nlm.ncbi.pubchemAxis2.TestedConcentrationType localTestedConcentration ;
                                
                           /*  This tracker boolean wil be used to detect whether the user called the set method
                          *   for this attribute. It will be used to determine whether to include this field
                           *   in the serialized XML
                           */
                           protected boolean localTestedConcentrationTracker = false ;
                           

                           /**
                           * Auto generated getter method
                           * @return gov.nih.nlm.ncbi.pubchemAxis2.TestedConcentrationType
                           */
                           public  gov.nih.nlm.ncbi.pubchemAxis2.TestedConcentrationType getTestedConcentration(){
                               return localTestedConcentration;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param TestedConcentration
                               */
                               public void setTestedConcentration(gov.nih.nlm.ncbi.pubchemAxis2.TestedConcentrationType param){
                            
                                       if (param != null){
                                          //update the setting tracker
                                          localTestedConcentrationTracker = true;
                                       } else {
                                          localTestedConcentrationTracker = false;
                                              
                                       }
                                   
                                            this.localTestedConcentration=param;
                                    

                               }
                            

                        /**
                        * field for ActiveConcentration
                        */

                        
                                    protected boolean localActiveConcentration ;
                                
                           /*  This tracker boolean wil be used to detect whether the user called the set method
                          *   for this attribute. It will be used to determine whether to include this field
                           *   in the serialized XML
                           */
                           protected boolean localActiveConcentrationTracker = false ;
                           

                           /**
                           * Auto generated getter method
                           * @return boolean
                           */
                           public  boolean getActiveConcentration(){
                               return localActiveConcentration;
                           }

                           
                        
                            /**
                               * Auto generated setter method
                               * @param param ActiveConcentration
                               */
                               public void setActiveConcentration(boolean param){
                            
                                       // setting primitive attribute tracker to true
                                       
                                               if (false) {
                                           localActiveConcentrationTracker = false;
                                              
                                       } else {
                                          localActiveConcentrationTracker = true;
                                       }
                                   
                                            this.localActiveConcentration=param;
                                    

                               }
                            

     /**
     * isReaderMTOMAware
     * @return true if the reader supports MTOM
     */
   public static boolean isReaderMTOMAware(javax.xml.stream.XMLStreamReader reader) {
        boolean isReaderMTOMAware = false;
        
        try{
          isReaderMTOMAware = java.lang.Boolean.TRUE.equals(reader.getProperty(org.apache.axiom.om.OMConstants.IS_DATA_HANDLERS_AWARE));
        }catch(java.lang.IllegalArgumentException e){
          isReaderMTOMAware = false;
        }
        return isReaderMTOMAware;
   }
     
     
        /**
        *
        * @param parentQName
        * @param factory
        * @return org.apache.axiom.om.OMElement
        */
       public org.apache.axiom.om.OMElement getOMElement (
               final javax.xml.namespace.QName parentQName,
               final org.apache.axiom.om.OMFactory factory) throws org.apache.axis2.databinding.ADBException{


        
               org.apache.axiom.om.OMDataSource dataSource =
                       new org.apache.axis2.databinding.ADBDataSource(this,parentQName){

                 public void serialize(org.apache.axis2.databinding.utils.writer.MTOMAwareXMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException {
                       ColumnDescriptionType.this.serialize(parentQName,factory,xmlWriter);
                 }
               };
               return new org.apache.axiom.om.impl.llom.OMSourcedElementImpl(
               parentQName,factory,dataSource);
            
       }

         public void serialize(final javax.xml.namespace.QName parentQName,
                                       final org.apache.axiom.om.OMFactory factory,
                                       org.apache.axis2.databinding.utils.writer.MTOMAwareXMLStreamWriter xmlWriter)
                                throws javax.xml.stream.XMLStreamException, org.apache.axis2.databinding.ADBException{
                           serialize(parentQName,factory,xmlWriter,false);
         }

         public void serialize(final javax.xml.namespace.QName parentQName,
                               final org.apache.axiom.om.OMFactory factory,
                               org.apache.axis2.databinding.utils.writer.MTOMAwareXMLStreamWriter xmlWriter,
                               boolean serializeType)
            throws javax.xml.stream.XMLStreamException, org.apache.axis2.databinding.ADBException{
            
                


                java.lang.String prefix = null;
                java.lang.String namespace = null;
                

                    prefix = parentQName.getPrefix();
                    namespace = parentQName.getNamespaceURI();

                    if ((namespace != null) && (namespace.trim().length() > 0)) {
                        java.lang.String writerPrefix = xmlWriter.getPrefix(namespace);
                        if (writerPrefix != null) {
                            xmlWriter.writeStartElement(namespace, parentQName.getLocalPart());
                        } else {
                            if (prefix == null) {
                                prefix = generatePrefix(namespace);
                            }

                            xmlWriter.writeStartElement(prefix, parentQName.getLocalPart(), namespace);
                            xmlWriter.writeNamespace(prefix, namespace);
                            xmlWriter.setPrefix(prefix, namespace);
                        }
                    } else {
                        xmlWriter.writeStartElement(parentQName.getLocalPart());
                    }
                
                  if (serializeType){
               

                   java.lang.String namespacePrefix = registerPrefix(xmlWriter,"http://pubchem.ncbi.nlm.nih.gov/");
                   if ((namespacePrefix != null) && (namespacePrefix.trim().length() > 0)){
                       writeAttribute("xsi","http://www.w3.org/2001/XMLSchema-instance","type",
                           namespacePrefix+":ColumnDescriptionType",
                           xmlWriter);
                   } else {
                       writeAttribute("xsi","http://www.w3.org/2001/XMLSchema-instance","type",
                           "ColumnDescriptionType",
                           xmlWriter);
                   }

               
                   }
               
                                            if (localHeading==null){
                                                 throw new org.apache.axis2.databinding.ADBException("Heading cannot be null!!");
                                            }
                                           localHeading.serialize(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Heading"),
                                               factory,xmlWriter);
                                         if (localTIDTracker){
                                    namespace = "http://pubchem.ncbi.nlm.nih.gov/";
                                    if (! namespace.equals("")) {
                                        prefix = xmlWriter.getPrefix(namespace);

                                        if (prefix == null) {
                                            prefix = generatePrefix(namespace);

                                            xmlWriter.writeStartElement(prefix,"TID", namespace);
                                            xmlWriter.writeNamespace(prefix, namespace);
                                            xmlWriter.setPrefix(prefix, namespace);

                                        } else {
                                            xmlWriter.writeStartElement(namespace,"TID");
                                        }

                                    } else {
                                        xmlWriter.writeStartElement("TID");
                                    }
                                
                                               if (localTID==java.lang.Integer.MIN_VALUE) {
                                           
                                                         throw new org.apache.axis2.databinding.ADBException("TID cannot be null!!");
                                                      
                                               } else {
                                                    xmlWriter.writeCharacters(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localTID));
                                               }
                                    
                                   xmlWriter.writeEndElement();
                             }
                                    namespace = "http://pubchem.ncbi.nlm.nih.gov/";
                                    if (! namespace.equals("")) {
                                        prefix = xmlWriter.getPrefix(namespace);

                                        if (prefix == null) {
                                            prefix = generatePrefix(namespace);

                                            xmlWriter.writeStartElement(prefix,"Name", namespace);
                                            xmlWriter.writeNamespace(prefix, namespace);
                                            xmlWriter.setPrefix(prefix, namespace);

                                        } else {
                                            xmlWriter.writeStartElement(namespace,"Name");
                                        }

                                    } else {
                                        xmlWriter.writeStartElement("Name");
                                    }
                                

                                          if (localName==null){
                                              // write the nil attribute
                                              
                                                     throw new org.apache.axis2.databinding.ADBException("Name cannot be null!!");
                                                  
                                          }else{

                                        
                                                   xmlWriter.writeCharacters(localName);
                                            
                                          }
                                    
                                   xmlWriter.writeEndElement();
                              if (localDescriptionTracker){
                                            if (localDescription==null){
                                                 throw new org.apache.axis2.databinding.ADBException("Description cannot be null!!");
                                            }
                                           localDescription.serialize(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Description"),
                                               factory,xmlWriter);
                                        }
                                    namespace = "http://pubchem.ncbi.nlm.nih.gov/";
                                    if (! namespace.equals("")) {
                                        prefix = xmlWriter.getPrefix(namespace);

                                        if (prefix == null) {
                                            prefix = generatePrefix(namespace);

                                            xmlWriter.writeStartElement(prefix,"Type", namespace);
                                            xmlWriter.writeNamespace(prefix, namespace);
                                            xmlWriter.setPrefix(prefix, namespace);

                                        } else {
                                            xmlWriter.writeStartElement(namespace,"Type");
                                        }

                                    } else {
                                        xmlWriter.writeStartElement("Type");
                                    }
                                

                                          if (localType==null){
                                              // write the nil attribute
                                              
                                                     throw new org.apache.axis2.databinding.ADBException("Type cannot be null!!");
                                                  
                                          }else{

                                        
                                                   xmlWriter.writeCharacters(localType);
                                            
                                          }
                                    
                                   xmlWriter.writeEndElement();
                              if (localUnitTracker){
                                    namespace = "http://pubchem.ncbi.nlm.nih.gov/";
                                    if (! namespace.equals("")) {
                                        prefix = xmlWriter.getPrefix(namespace);

                                        if (prefix == null) {
                                            prefix = generatePrefix(namespace);

                                            xmlWriter.writeStartElement(prefix,"Unit", namespace);
                                            xmlWriter.writeNamespace(prefix, namespace);
                                            xmlWriter.setPrefix(prefix, namespace);

                                        } else {
                                            xmlWriter.writeStartElement(namespace,"Unit");
                                        }

                                    } else {
                                        xmlWriter.writeStartElement("Unit");
                                    }
                                

                                          if (localUnit==null){
                                              // write the nil attribute
                                              
                                                     throw new org.apache.axis2.databinding.ADBException("Unit cannot be null!!");
                                                  
                                          }else{

                                        
                                                   xmlWriter.writeCharacters(localUnit);
                                            
                                          }
                                    
                                   xmlWriter.writeEndElement();
                             } if (localTestedConcentrationTracker){
                                            if (localTestedConcentration==null){
                                                 throw new org.apache.axis2.databinding.ADBException("TestedConcentration cannot be null!!");
                                            }
                                           localTestedConcentration.serialize(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","TestedConcentration"),
                                               factory,xmlWriter);
                                        } if (localActiveConcentrationTracker){
                                    namespace = "http://pubchem.ncbi.nlm.nih.gov/";
                                    if (! namespace.equals("")) {
                                        prefix = xmlWriter.getPrefix(namespace);

                                        if (prefix == null) {
                                            prefix = generatePrefix(namespace);

                                            xmlWriter.writeStartElement(prefix,"ActiveConcentration", namespace);
                                            xmlWriter.writeNamespace(prefix, namespace);
                                            xmlWriter.setPrefix(prefix, namespace);

                                        } else {
                                            xmlWriter.writeStartElement(namespace,"ActiveConcentration");
                                        }

                                    } else {
                                        xmlWriter.writeStartElement("ActiveConcentration");
                                    }
                                
                                               if (false) {
                                           
                                                         throw new org.apache.axis2.databinding.ADBException("ActiveConcentration cannot be null!!");
                                                      
                                               } else {
                                                    xmlWriter.writeCharacters(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localActiveConcentration));
                                               }
                                    
                                   xmlWriter.writeEndElement();
                             }
                    xmlWriter.writeEndElement();
               

        }

         /**
          * Util method to write an attribute with the ns prefix
          */
          private void writeAttribute(java.lang.String prefix,java.lang.String namespace,java.lang.String attName,
                                      java.lang.String attValue,javax.xml.stream.XMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException{
              if (xmlWriter.getPrefix(namespace) == null) {
                       xmlWriter.writeNamespace(prefix, namespace);
                       xmlWriter.setPrefix(prefix, namespace);

              }

              xmlWriter.writeAttribute(namespace,attName,attValue);

         }

        /**
          * Util method to write an attribute without the ns prefix
          */
          private void writeAttribute(java.lang.String namespace,java.lang.String attName,
                                      java.lang.String attValue,javax.xml.stream.XMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException{
                if (namespace.equals(""))
              {
                  xmlWriter.writeAttribute(attName,attValue);
              }
              else
              {
                  registerPrefix(xmlWriter, namespace);
                  xmlWriter.writeAttribute(namespace,attName,attValue);
              }
          }


           /**
             * Util method to write an attribute without the ns prefix
             */
            private void writeQNameAttribute(java.lang.String namespace, java.lang.String attName,
                                             javax.xml.namespace.QName qname, javax.xml.stream.XMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException {

                java.lang.String attributeNamespace = qname.getNamespaceURI();
                java.lang.String attributePrefix = xmlWriter.getPrefix(attributeNamespace);
                if (attributePrefix == null) {
                    attributePrefix = registerPrefix(xmlWriter, attributeNamespace);
                }
                java.lang.String attributeValue;
                if (attributePrefix.trim().length() > 0) {
                    attributeValue = attributePrefix + ":" + qname.getLocalPart();
                } else {
                    attributeValue = qname.getLocalPart();
                }

                if (namespace.equals("")) {
                    xmlWriter.writeAttribute(attName, attributeValue);
                } else {
                    registerPrefix(xmlWriter, namespace);
                    xmlWriter.writeAttribute(namespace, attName, attributeValue);
                }
            }
        /**
         *  method to handle Qnames
         */

        private void writeQName(javax.xml.namespace.QName qname,
                                javax.xml.stream.XMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException {
            java.lang.String namespaceURI = qname.getNamespaceURI();
            if (namespaceURI != null) {
                java.lang.String prefix = xmlWriter.getPrefix(namespaceURI);
                if (prefix == null) {
                    prefix = generatePrefix(namespaceURI);
                    xmlWriter.writeNamespace(prefix, namespaceURI);
                    xmlWriter.setPrefix(prefix,namespaceURI);
                }

                if (prefix.trim().length() > 0){
                    xmlWriter.writeCharacters(prefix + ":" + org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qname));
                } else {
                    // i.e this is the default namespace
                    xmlWriter.writeCharacters(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qname));
                }

            } else {
                xmlWriter.writeCharacters(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qname));
            }
        }

        private void writeQNames(javax.xml.namespace.QName[] qnames,
                                 javax.xml.stream.XMLStreamWriter xmlWriter) throws javax.xml.stream.XMLStreamException {

            if (qnames != null) {
                // we have to store this data until last moment since it is not possible to write any
                // namespace data after writing the charactor data
                java.lang.StringBuffer stringToWrite = new java.lang.StringBuffer();
                java.lang.String namespaceURI = null;
                java.lang.String prefix = null;

                for (int i = 0; i < qnames.length; i++) {
                    if (i > 0) {
                        stringToWrite.append(" ");
                    }
                    namespaceURI = qnames[i].getNamespaceURI();
                    if (namespaceURI != null) {
                        prefix = xmlWriter.getPrefix(namespaceURI);
                        if ((prefix == null) || (prefix.length() == 0)) {
                            prefix = generatePrefix(namespaceURI);
                            xmlWriter.writeNamespace(prefix, namespaceURI);
                            xmlWriter.setPrefix(prefix,namespaceURI);
                        }

                        if (prefix.trim().length() > 0){
                            stringToWrite.append(prefix).append(":").append(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qnames[i]));
                        } else {
                            stringToWrite.append(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qnames[i]));
                        }
                    } else {
                        stringToWrite.append(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(qnames[i]));
                    }
                }
                xmlWriter.writeCharacters(stringToWrite.toString());
            }

        }


         /**
         * Register a namespace prefix
         */
         private java.lang.String registerPrefix(javax.xml.stream.XMLStreamWriter xmlWriter, java.lang.String namespace) throws javax.xml.stream.XMLStreamException {
                java.lang.String prefix = xmlWriter.getPrefix(namespace);

                if (prefix == null) {
                    prefix = generatePrefix(namespace);

                    while (xmlWriter.getNamespaceContext().getNamespaceURI(prefix) != null) {
                        prefix = org.apache.axis2.databinding.utils.BeanUtil.getUniquePrefix();
                    }

                    xmlWriter.writeNamespace(prefix, namespace);
                    xmlWriter.setPrefix(prefix, namespace);
                }

                return prefix;
            }


  
        /**
        * databinding method to get an XML representation of this object
        *
        */
        public javax.xml.stream.XMLStreamReader getPullParser(javax.xml.namespace.QName qName)
                    throws org.apache.axis2.databinding.ADBException{


        
                 java.util.ArrayList elementList = new java.util.ArrayList();
                 java.util.ArrayList attribList = new java.util.ArrayList();

                
                            elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "Heading"));
                            
                            
                                    if (localHeading==null){
                                         throw new org.apache.axis2.databinding.ADBException("Heading cannot be null!!");
                                    }
                                    elementList.add(localHeading);
                                 if (localTIDTracker){
                                      elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "TID"));
                                 
                                elementList.add(
                                   org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localTID));
                            }
                                      elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "Name"));
                                 
                                        if (localName != null){
                                            elementList.add(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localName));
                                        } else {
                                           throw new org.apache.axis2.databinding.ADBException("Name cannot be null!!");
                                        }
                                     if (localDescriptionTracker){
                            elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "Description"));
                            
                            
                                    if (localDescription==null){
                                         throw new org.apache.axis2.databinding.ADBException("Description cannot be null!!");
                                    }
                                    elementList.add(localDescription);
                                }
                                      elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "Type"));
                                 
                                        if (localType != null){
                                            elementList.add(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localType));
                                        } else {
                                           throw new org.apache.axis2.databinding.ADBException("Type cannot be null!!");
                                        }
                                     if (localUnitTracker){
                                      elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "Unit"));
                                 
                                        if (localUnit != null){
                                            elementList.add(org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localUnit));
                                        } else {
                                           throw new org.apache.axis2.databinding.ADBException("Unit cannot be null!!");
                                        }
                                    } if (localTestedConcentrationTracker){
                            elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "TestedConcentration"));
                            
                            
                                    if (localTestedConcentration==null){
                                         throw new org.apache.axis2.databinding.ADBException("TestedConcentration cannot be null!!");
                                    }
                                    elementList.add(localTestedConcentration);
                                } if (localActiveConcentrationTracker){
                                      elementList.add(new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/",
                                                                      "ActiveConcentration"));
                                 
                                elementList.add(
                                   org.apache.axis2.databinding.utils.ConverterUtil.convertToString(localActiveConcentration));
                            }

                return new org.apache.axis2.databinding.utils.reader.ADBXMLStreamReaderImpl(qName, elementList.toArray(), attribList.toArray());
            
            

        }

  

     /**
      *  Factory class that keeps the parse method
      */
    public static class Factory{

        
        

        /**
        * static method to create the object
        * Precondition:  If this object is an element, the current or next start element starts this object and any intervening reader events are ignorable
        *                If this object is not an element, it is a complex type and the reader is at the event just after the outer start element
        * Postcondition: If this object is an element, the reader is positioned at its end element
        *                If this object is a complex type, the reader is positioned at the end element of its outer element
        */
        public static ColumnDescriptionType parse(javax.xml.stream.XMLStreamReader reader) throws java.lang.Exception{
            ColumnDescriptionType object =
                new ColumnDescriptionType();

            int event;
            java.lang.String nillableValue = null;
            java.lang.String prefix ="";
            java.lang.String namespaceuri ="";
            try {
                
                while (!reader.isStartElement() && !reader.isEndElement())
                    reader.next();

                
                if (reader.getAttributeValue("http://www.w3.org/2001/XMLSchema-instance","type")!=null){
                  java.lang.String fullTypeName = reader.getAttributeValue("http://www.w3.org/2001/XMLSchema-instance",
                        "type");
                  if (fullTypeName!=null){
                    java.lang.String nsPrefix = null;
                    if (fullTypeName.indexOf(":") > -1){
                        nsPrefix = fullTypeName.substring(0,fullTypeName.indexOf(":"));
                    }
                    nsPrefix = nsPrefix==null?"":nsPrefix;

                    java.lang.String type = fullTypeName.substring(fullTypeName.indexOf(":")+1);
                    
                            if (!"ColumnDescriptionType".equals(type)){
                                //find namespace for the prefix
                                java.lang.String nsUri = reader.getNamespaceContext().getNamespaceURI(nsPrefix);
                                return (ColumnDescriptionType)gov.nih.nlm.ncbi.pubchemAxis2.ExtensionMapper.getTypeObject(
                                     nsUri,type,reader);
                              }
                        

                  }
                

                }

                

                
                // Note all attributes that were handled. Used to differ normal attributes
                // from anyAttributes.
                java.util.Vector handledAttributes = new java.util.Vector();
                

                 
                    
                    reader.next();
                
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Heading").equals(reader.getName())){
                                
                                                object.setHeading(gov.nih.nlm.ncbi.pubchemAxis2.HeadingType.Factory.parse(reader));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                else{
                                    // A start element we are not expecting indicates an invalid parameter was passed
                                    throw new org.apache.axis2.databinding.ADBException("Unexpected subelement " + reader.getLocalName());
                                }
                            
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","TID").equals(reader.getName())){
                                
                                    java.lang.String content = reader.getElementText();
                                    
                                              object.setTID(
                                                    org.apache.axis2.databinding.utils.ConverterUtil.convertToInt(content));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                    else {
                                        
                                               object.setTID(java.lang.Integer.MIN_VALUE);
                                           
                                    }
                                
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Name").equals(reader.getName())){
                                
                                    java.lang.String content = reader.getElementText();
                                    
                                              object.setName(
                                                    org.apache.axis2.databinding.utils.ConverterUtil.convertToString(content));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                else{
                                    // A start element we are not expecting indicates an invalid parameter was passed
                                    throw new org.apache.axis2.databinding.ADBException("Unexpected subelement " + reader.getLocalName());
                                }
                            
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Description").equals(reader.getName())){
                                
                                                object.setDescription(gov.nih.nlm.ncbi.pubchemAxis2.ArrayOfString.Factory.parse(reader));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                    else {
                                        
                                    }
                                
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Type").equals(reader.getName())){
                                
                                    java.lang.String content = reader.getElementText();
                                    
                                              object.setType(
                                                    org.apache.axis2.databinding.utils.ConverterUtil.convertToString(content));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                else{
                                    // A start element we are not expecting indicates an invalid parameter was passed
                                    throw new org.apache.axis2.databinding.ADBException("Unexpected subelement " + reader.getLocalName());
                                }
                            
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","Unit").equals(reader.getName())){
                                
                                    java.lang.String content = reader.getElementText();
                                    
                                              object.setUnit(
                                                    org.apache.axis2.databinding.utils.ConverterUtil.convertToString(content));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                    else {
                                        
                                    }
                                
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","TestedConcentration").equals(reader.getName())){
                                
                                                object.setTestedConcentration(gov.nih.nlm.ncbi.pubchemAxis2.TestedConcentrationType.Factory.parse(reader));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                    else {
                                        
                                    }
                                
                                    
                                    while (!reader.isStartElement() && !reader.isEndElement()) reader.next();
                                
                                    if (reader.isStartElement() && new javax.xml.namespace.QName("http://pubchem.ncbi.nlm.nih.gov/","ActiveConcentration").equals(reader.getName())){
                                
                                    java.lang.String content = reader.getElementText();
                                    
                                              object.setActiveConcentration(
                                                    org.apache.axis2.databinding.utils.ConverterUtil.convertToBoolean(content));
                                              
                                        reader.next();
                                    
                              }  // End of if for expected property start element
                                
                                    else {
                                        
                                    }
                                  
                            while (!reader.isStartElement() && !reader.isEndElement())
                                reader.next();
                            
                                if (reader.isStartElement())
                                // A start element we are not expecting indicates a trailing invalid property
                                throw new org.apache.axis2.databinding.ADBException("Unexpected subelement " + reader.getLocalName());
                            



            } catch (javax.xml.stream.XMLStreamException e) {
                throw new java.lang.Exception(e);
            }

            return object;
        }

        }//end of factory class

        

        }
           
          