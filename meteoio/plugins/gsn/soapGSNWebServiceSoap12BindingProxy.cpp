/* soapGSNWebServiceSoap12BindingProxy.cpp
   Generated by gSOAP 2.7.9l from GSNWebService.h
   Copyright(C) 2000-2007, Robert van Engelen, Genivia Inc. All Rights Reserved.
   This part of the software is released under one of the following licenses:
   GPL, the gSOAP public license, or Genivia's license for commercial use.
*/

#include "soapGSNWebServiceSoap12BindingProxy.h"

GSNWebServiceSoap12BindingProxy::GSNWebServiceSoap12BindingProxy()
{	GSNWebServiceSoap12BindingProxy_init(SOAP_IO_DEFAULT, SOAP_IO_DEFAULT);
}

GSNWebServiceSoap12BindingProxy::GSNWebServiceSoap12BindingProxy(soap_mode iomode)
{	GSNWebServiceSoap12BindingProxy_init(iomode, iomode);
}

GSNWebServiceSoap12BindingProxy::GSNWebServiceSoap12BindingProxy(soap_mode imode, soap_mode omode)
{	GSNWebServiceSoap12BindingProxy_init(imode, omode);
}

void GSNWebServiceSoap12BindingProxy::GSNWebServiceSoap12BindingProxy_init(soap_mode imode, soap_mode omode)
{	soap_imode(this, imode);
	soap_omode(this, omode);
	soap_endpoint = NULL;
	static const struct Namespace namespaces[] =
{
	{"SOAP-ENV", "http://www.w3.org/2003/05/soap-envelope", "http://www.w3.org/*/soap-envelope", NULL},
	{"SOAP-ENC", "http://www.w3.org/2003/05/soap-encoding", "http://www.w3.org/*/soap-encoding", NULL},
	{"xsi", "http://www.w3.org/2001/XMLSchema-instance", "http://www.w3.org/*/XMLSchema-instance", NULL},
	{"xsd", "http://www.w3.org/2001/XMLSchema", "http://www.w3.org/*/XMLSchema", NULL},
	{"ns3", "http://datarequest.http.gsn/xsd", NULL, NULL},
	{"ns2", "http://standard.webservice.gsn/xsd", NULL, NULL},
	{"ns4", "http://standard.webservice.gsn/GSNWebServiceSoap12Binding", NULL, NULL},
	{"ns1", "http://standard.webservice.gsn", NULL, NULL},
	{"ns5", "http://standard.webservice.gsn/GSNWebServiceSoap11Binding", NULL, NULL},
	{NULL, NULL, NULL, NULL}
};
	if (!this->namespaces)
		this->namespaces = namespaces;
}

GSNWebServiceSoap12BindingProxy::~GSNWebServiceSoap12BindingProxy()
{ }

void GSNWebServiceSoap12BindingProxy::soap_noheader()
{	header = NULL;
}

const SOAP_ENV__Fault *GSNWebServiceSoap12BindingProxy::soap_fault()
{	return this->fault;
}

const char *GSNWebServiceSoap12BindingProxy::soap_fault_string()
{	return *soap_faultstring(this);
}

const char *GSNWebServiceSoap12BindingProxy::soap_fault_detail()
{	return *soap_faultdetail(this);
}

int GSNWebServiceSoap12BindingProxy::getContainerInfo(_ns1__getContainerInfo *ns1__getContainerInfo, _ns1__getContainerInfoResponse *ns1__getContainerInfoResponse)
{	struct soap *soap = this;
	struct __ns4__getContainerInfo soap_tmp___ns4__getContainerInfo;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:getContainerInfo";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__getContainerInfo.ns1__getContainerInfo = ns1__getContainerInfo;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__getContainerInfo(soap, &soap_tmp___ns4__getContainerInfo);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__getContainerInfo(soap, &soap_tmp___ns4__getContainerInfo, "-ns4:getContainerInfo", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__getContainerInfo(soap, &soap_tmp___ns4__getContainerInfo, "-ns4:getContainerInfo", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__getContainerInfoResponse)
		return soap_closesock(soap);
	ns1__getContainerInfoResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__getContainerInfoResponse->soap_get(soap, "ns1:getContainerInfoResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::getMultiData(_ns1__getMultiData *ns1__getMultiData, _ns1__getMultiDataResponse *ns1__getMultiDataResponse)
{	struct soap *soap = this;
	struct __ns4__getMultiData soap_tmp___ns4__getMultiData;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:getMultiData";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__getMultiData.ns1__getMultiData = ns1__getMultiData;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__getMultiData(soap, &soap_tmp___ns4__getMultiData);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__getMultiData(soap, &soap_tmp___ns4__getMultiData, "-ns4:getMultiData", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__getMultiData(soap, &soap_tmp___ns4__getMultiData, "-ns4:getMultiData", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__getMultiDataResponse)
		return soap_closesock(soap);
	ns1__getMultiDataResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__getMultiDataResponse->soap_get(soap, "ns1:getMultiDataResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::listVirtualSensorNames(_ns1__listVirtualSensorNames *ns1__listVirtualSensorNames, _ns1__listVirtualSensorNamesResponse *ns1__listVirtualSensorNamesResponse)
{	struct soap *soap = this;
	struct __ns4__listVirtualSensorNames soap_tmp___ns4__listVirtualSensorNames;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:listVirtualSensorNames";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__listVirtualSensorNames.ns1__listVirtualSensorNames = ns1__listVirtualSensorNames;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__listVirtualSensorNames(soap, &soap_tmp___ns4__listVirtualSensorNames);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__listVirtualSensorNames(soap, &soap_tmp___ns4__listVirtualSensorNames, "-ns4:listVirtualSensorNames", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__listVirtualSensorNames(soap, &soap_tmp___ns4__listVirtualSensorNames, "-ns4:listVirtualSensorNames", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__listVirtualSensorNamesResponse)
		return soap_closesock(soap);
	ns1__listVirtualSensorNamesResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__listVirtualSensorNamesResponse->soap_get(soap, "ns1:listVirtualSensorNamesResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::getNextData(_ns1__getNextData *ns1__getNextData, _ns1__getNextDataResponse *ns1__getNextDataResponse)
{	struct soap *soap = this;
	struct __ns4__getNextData soap_tmp___ns4__getNextData;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:getNextData";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__getNextData.ns1__getNextData = ns1__getNextData;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__getNextData(soap, &soap_tmp___ns4__getNextData);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__getNextData(soap, &soap_tmp___ns4__getNextData, "-ns4:getNextData", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__getNextData(soap, &soap_tmp___ns4__getNextData, "-ns4:getNextData", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__getNextDataResponse)
		return soap_closesock(soap);
	ns1__getNextDataResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__getNextDataResponse->soap_get(soap, "ns1:getNextDataResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::registerQuery(_ns1__registerQuery *ns1__registerQuery, _ns1__registerQueryResponse *ns1__registerQueryResponse)
{	struct soap *soap = this;
	struct __ns4__registerQuery soap_tmp___ns4__registerQuery;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:registerQuery";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__registerQuery.ns1__registerQuery = ns1__registerQuery;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__registerQuery(soap, &soap_tmp___ns4__registerQuery);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__registerQuery(soap, &soap_tmp___ns4__registerQuery, "-ns4:registerQuery", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__registerQuery(soap, &soap_tmp___ns4__registerQuery, "-ns4:registerQuery", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__registerQueryResponse)
		return soap_closesock(soap);
	ns1__registerQueryResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__registerQueryResponse->soap_get(soap, "ns1:registerQueryResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::unregisterQuery(_ns1__unregisterQuery *ns1__unregisterQuery, _ns1__unregisterQueryResponse *ns1__unregisterQueryResponse)
{	struct soap *soap = this;
	struct __ns4__unregisterQuery soap_tmp___ns4__unregisterQuery;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:unregisterQuery";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__unregisterQuery.ns1__unregisterQuery = ns1__unregisterQuery;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__unregisterQuery(soap, &soap_tmp___ns4__unregisterQuery);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__unregisterQuery(soap, &soap_tmp___ns4__unregisterQuery, "-ns4:unregisterQuery", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__unregisterQuery(soap, &soap_tmp___ns4__unregisterQuery, "-ns4:unregisterQuery", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__unregisterQueryResponse)
		return soap_closesock(soap);
	ns1__unregisterQueryResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__unregisterQueryResponse->soap_get(soap, "ns1:unregisterQueryResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::getLatestMultiData(_ns1__getLatestMultiData *ns1__getLatestMultiData, _ns1__getLatestMultiDataResponse *ns1__getLatestMultiDataResponse)
{	struct soap *soap = this;
	struct __ns4__getLatestMultiData soap_tmp___ns4__getLatestMultiData;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:getLatestMultiData";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__getLatestMultiData.ns1__getLatestMultiData = ns1__getLatestMultiData;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__getLatestMultiData(soap, &soap_tmp___ns4__getLatestMultiData);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__getLatestMultiData(soap, &soap_tmp___ns4__getLatestMultiData, "-ns4:getLatestMultiData", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__getLatestMultiData(soap, &soap_tmp___ns4__getLatestMultiData, "-ns4:getLatestMultiData", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__getLatestMultiDataResponse)
		return soap_closesock(soap);
	ns1__getLatestMultiDataResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__getLatestMultiDataResponse->soap_get(soap, "ns1:getLatestMultiDataResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::deleteVirtualSensor(_ns1__deleteVirtualSensor *ns1__deleteVirtualSensor, _ns1__deleteVirtualSensorResponse *ns1__deleteVirtualSensorResponse)
{	struct soap *soap = this;
	struct __ns4__deleteVirtualSensor soap_tmp___ns4__deleteVirtualSensor;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:deleteVirtualSensor";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__deleteVirtualSensor.ns1__deleteVirtualSensor = ns1__deleteVirtualSensor;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__deleteVirtualSensor(soap, &soap_tmp___ns4__deleteVirtualSensor);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__deleteVirtualSensor(soap, &soap_tmp___ns4__deleteVirtualSensor, "-ns4:deleteVirtualSensor", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__deleteVirtualSensor(soap, &soap_tmp___ns4__deleteVirtualSensor, "-ns4:deleteVirtualSensor", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__deleteVirtualSensorResponse)
		return soap_closesock(soap);
	ns1__deleteVirtualSensorResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__deleteVirtualSensorResponse->soap_get(soap, "ns1:deleteVirtualSensorResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::getVirtualSensorsDetails(_ns1__getVirtualSensorsDetails *ns1__getVirtualSensorsDetails, _ns1__getVirtualSensorsDetailsResponse *ns1__getVirtualSensorsDetailsResponse)
{	struct soap *soap = this;
	struct __ns4__getVirtualSensorsDetails soap_tmp___ns4__getVirtualSensorsDetails;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:getVirtualSensorsDetails";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__getVirtualSensorsDetails.ns1__getVirtualSensorsDetails = ns1__getVirtualSensorsDetails;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__getVirtualSensorsDetails(soap, &soap_tmp___ns4__getVirtualSensorsDetails);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__getVirtualSensorsDetails(soap, &soap_tmp___ns4__getVirtualSensorsDetails, "-ns4:getVirtualSensorsDetails", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__getVirtualSensorsDetails(soap, &soap_tmp___ns4__getVirtualSensorsDetails, "-ns4:getVirtualSensorsDetails", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__getVirtualSensorsDetailsResponse)
		return soap_closesock(soap);
	ns1__getVirtualSensorsDetailsResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__getVirtualSensorsDetailsResponse->soap_get(soap, "ns1:getVirtualSensorsDetailsResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::createVirtualSensor(_ns1__createVirtualSensor *ns1__createVirtualSensor, _ns1__createVirtualSensorResponse *ns1__createVirtualSensorResponse)
{	struct soap *soap = this;
	struct __ns4__createVirtualSensor soap_tmp___ns4__createVirtualSensor;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:createVirtualSensor";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__createVirtualSensor.ns1__createVirtualSensor = ns1__createVirtualSensor;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__createVirtualSensor(soap, &soap_tmp___ns4__createVirtualSensor);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__createVirtualSensor(soap, &soap_tmp___ns4__createVirtualSensor, "-ns4:createVirtualSensor", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__createVirtualSensor(soap, &soap_tmp___ns4__createVirtualSensor, "-ns4:createVirtualSensor", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__createVirtualSensorResponse)
		return soap_closesock(soap);
	ns1__createVirtualSensorResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__createVirtualSensorResponse->soap_get(soap, "ns1:createVirtualSensorResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}

int GSNWebServiceSoap12BindingProxy::listWrapperURLs(_ns1__listWrapperURLs *ns1__listWrapperURLs, _ns1__listWrapperURLsResponse *ns1__listWrapperURLsResponse)
{	struct soap *soap = this;
	struct __ns4__listWrapperURLs soap_tmp___ns4__listWrapperURLs;
	const char *soap_action = NULL;
	if (!soap_endpoint)
		soap_endpoint = "http://192.33.210.37:22001/services/GSNWebService/";
	soap_action = "urn:listWrapperURLs";
	soap->encodingStyle = NULL;
	soap_tmp___ns4__listWrapperURLs.ns1__listWrapperURLs = ns1__listWrapperURLs;
	soap_begin(soap);
	soap_serializeheader(soap);
	soap_serialize___ns4__listWrapperURLs(soap, &soap_tmp___ns4__listWrapperURLs);
	if (soap_begin_count(soap))
		return soap->error;
	if (soap->mode & SOAP_IO_LENGTH)
	{	if (soap_envelope_begin_out(soap)
		 || soap_putheader(soap)
		 || soap_body_begin_out(soap)
		 || soap_put___ns4__listWrapperURLs(soap, &soap_tmp___ns4__listWrapperURLs, "-ns4:listWrapperURLs", "")
		 || soap_body_end_out(soap)
		 || soap_envelope_end_out(soap))
			 return soap->error;
	}
	if (soap_end_count(soap))
		return soap->error;
	if (soap_connect(soap, soap_endpoint, soap_action)
	 || soap_envelope_begin_out(soap)
	 || soap_putheader(soap)
	 || soap_body_begin_out(soap)
	 || soap_put___ns4__listWrapperURLs(soap, &soap_tmp___ns4__listWrapperURLs, "-ns4:listWrapperURLs", "")
	 || soap_body_end_out(soap)
	 || soap_envelope_end_out(soap)
	 || soap_end_send(soap))
		return soap_closesock(soap);
	if (!ns1__listWrapperURLsResponse)
		return soap_closesock(soap);
	ns1__listWrapperURLsResponse->soap_default(soap);
	if (soap_begin_recv(soap)
	 || soap_envelope_begin_in(soap)
	 || soap_recv_header(soap)
	 || soap_body_begin_in(soap))
		return soap_closesock(soap);
	ns1__listWrapperURLsResponse->soap_get(soap, "ns1:listWrapperURLsResponse", "");
	if (soap->error)
	{	if (soap->error == SOAP_TAG_MISMATCH && soap->level == 2)
			return soap_recv_fault(soap);
		return soap_closesock(soap);
	}
	if (soap_body_end_in(soap)
	 || soap_envelope_end_in(soap)
	 || soap_end_recv(soap))
		return soap_closesock(soap);
	return soap_closesock(soap);
}
/* End of client proxy code */
