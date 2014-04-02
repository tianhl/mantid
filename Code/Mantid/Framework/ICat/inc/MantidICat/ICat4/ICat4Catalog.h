#ifndef MANTID_ICAT_ICAT4CATALOG_H_
#define MANTID_ICAT_ICAT4CATALOG_H_

#include "MantidAPI/ICatalog.h"
#include "MantidAPI/ICatalogInfoService.h"
#include "MantidAPI/TableRow.h"
#include "MantidICat/CatalogHelper.h"
#include "MantidICat/CatalogSearchParam.h"

namespace Mantid
{
  namespace ICat
  {
    /**
    This class is responsible for the implementation of ICat4 version based information catalogs
    @author Jay Rainey, ISIS Rutherford Appleton Laboratory 
    @date 12/08/2013
    Copyright &copy; 2013 ISIS Rutherford Appleton Laboratory & NScD Oak Ridge National Laboratory

    This file is part of Mantid.

    Mantid is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Mantid is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    File change history is stored at: <https://github.com/mantidproject/mantid>.
    Code Documentation is available at: <http://doxygen.mantidproject.org>
    */
    class ICat4Catalog : public Mantid::API::ICatalog, public Mantid::API::ICatalogInfoService
    {
      public:
        /// Constructor
        ICat4Catalog();

        /// Log the user into the catalog system.
        virtual API::CatalogSession_sptr login(const std::string& username,const std::string& password,
            const std::string& endpoint,const std::string& facility);
        /// Log the user out of the catalog system.
        virtual void logout();
        /// Search the catalog for data.
        virtual void search(const CatalogSearchParam& inputs, Mantid::API::ITableWorkspace_sptr& outputws,
            const int &offset, const int &limit);
        /// Obtain the number of results returned by the search method.
        virtual int64_t getNumberOfSearchResults(const CatalogSearchParam& inputs);
        /// Show the logged in user's investigations search results.
        virtual void myData(Mantid::API::ITableWorkspace_sptr& outputws);
        /// Get datasets.
        virtual void getDataSets(const std::string&investigationId,Mantid::API::ITableWorkspace_sptr& outputws);
        /// Get datafiles
        virtual void getDataFiles(const std::string&investigationId,Mantid::API::ITableWorkspace_sptr& outputws);
        /// Get instruments list
        virtual void listInstruments(std::vector<std::string>& instruments);
        /// Get investigationtypes list
        virtual void listInvestigationTypes(std::vector<std::string>& invstTypes);
        /// Keep current session alive
        virtual void keepAlive();

        /// Get the file location string(s) from archive.
        virtual const std::string getFileLocation(const long long&fileID);
        /// Get the url(s) based on the fileID.
        virtual const std::string getDownloadURL(const long long& fileID);
        /// get URL of where to PUT (publish) files.
        virtual const std::string getUploadURL(
            const std::string &investigationID, const std::string &createFileName, const std::string &dataFileDescription);

      private:
        // Defines the SSL authentication scheme.
        void setSSLContext(ICat4::ICATPortBindingProxy& icat);
        // Creates a search query string based on inputs provided by the user.
        std::string buildSearchQuery(const CatalogSearchParam& inputs);
        // Search the archive & obtain the dataset ID based on the investigationID.
        int64_t getDatasetId(const std::string &investigationID);
        // Sets the soap-endpoint & SSL context for the given ICAT proxy.
        void setICATProxySettings(ICat4::ICATPortBindingProxy& icat);

        // Allows easy access to helper methods.
        CatalogHelper m_catalogHelper;
    };
  }
}

#endif
