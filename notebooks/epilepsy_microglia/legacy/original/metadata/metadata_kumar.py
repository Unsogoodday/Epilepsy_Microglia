from typing import Dict, TypedDict, Literal

DATASET = "GSE201048"
PROTOCOL = "CITE-seq"

Dx = Literal["OLE", "FCD", "TLE", "SWS", "SRFSE"]
Region = Literal["Frontal", "Temporal", "Occipital"]
Hemisphere = Literal["L", "R"]
Sex = Literal["M", "F"]


class SampleMetadata(TypedDict):
    patient_id: str
    dataset: str
    sex: Sex
    age: int
    dx: Dx
    dx_subtype: str
    region: Region
    hemisphere: Hemisphere
    procedure: str
    protocol: str


def kumar_metadata() -> Dict[str, SampleMetadata]:
    return {
        "GSM6049632": {
            "patient_id": "Kumar01A",
            "dataset": DATASET,
            "sex": "F",
            "age": 9,
            "dx": "OLE",
            "dx_subtype": "OLE",
            "region": "Occipital",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
        },
        "GSM6049633": {
            "patient_id": "Kumar01B",
            "dataset": DATASET,
            "sex": "F",
            "age": 9,
            "dx": "OLE",
            "dx_subtype": "OLE",
            "region": "Occipital",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
        },
        "GSM6049634": {
            "patient_id": "Kumar02",
            "dataset": DATASET,
            "sex": "F",
            "age": 4,
            "dx": "FCD",
            "dx_subtype": "FCD_IIb",
            "region": "Frontal",
            "hemisphere": "R",
            "procedure": "lesionectomy",
            "protocol": PROTOCOL,
        },
                "GSM6049635": {
            "patient_id": "Kumar03",
            "dataset": DATASET,
            "sex": "M",
            "age": 18,
            "dx": "FCD",
            "dx_subtype": "FCD_IIb",
            "region": "Frontal",
            "hemisphere": "L",
            "procedure": "lesionectomy",
            "protocol": PROTOCOL,
            "replicate": "A",
        },
        "GSM6049636": {
            "patient_id": "Kumar03",
            "dataset": DATASET,
            "sex": "M",
            "age": 18,
            "dx": "FCD",
            "dx_subtype": "FCD_IIb",
            "region": "Frontal",
            "hemisphere": "L",
            "procedure": "lesionectomy",
            "protocol": PROTOCOL,
            "replicate": "B",
        },
        "GSM6049637": {
            "patient_id": "Kumar04",
            "dataset": DATASET,
            "sex": "F",
            "age": 3,
            "dx": "TLE",
            "dx_subtype": "TLE",
            "region": "Temporal",
            "hemisphere": "L",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
        },
        "GSM6049638": {
            "patient_id": "Kumar05",
            "dataset": DATASET,
            "sex": "F",
            "age": 22,
            "dx": "SWS",
            "dx_subtype": "SWS",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "A",
        },
        "GSM6049639": {
            "patient_id": "Kumar05",
            "dataset": DATASET,
            "sex": "F",
            "age": 22,
            "dx": "SWS",
            "dx_subtype": "SWS",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "B",
        },
        "GSM6049640": {
            "patient_id": "Kumar06",
            "dataset": DATASET,
            "sex": "F",
            "age": 4,
            "dx": "SRFSE",
            "dx_subtype": "SRFSE",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "A",
        },
        "GSM6049641": {
            "patient_id": "Kumar06",
            "dataset": DATASET,
            "sex": "F",
            "age": 4,
            "dx": "SRFSE",
            "dx_subtype": "SRFSE",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "B1",
        },
        "GSM6049642": {
            "patient_id": "Kumar06",
            "dataset": DATASET,
            "sex": "F",
            "age": 4,
            "dx": "SRFSE",
            "dx_subtype": "SRFSE",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "B2",
        },
        "GSM6049643": {
            "patient_id": "Kumar06",
            "dataset": DATASET,
            "sex": "F",
            "age": 4,
            "dx": "SRFSE",
            "dx_subtype": "SRFSE",
            "region": "Temporal",
            "hemisphere": "R",
            "procedure": "lobectomy",
            "protocol": PROTOCOL,
            "replicate": "B3",
        },
    }
