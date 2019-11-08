from GUD import GUDUtils
from GUD.ORM import ClinVar
db_name = GUDUtils._get_db_name()
engine, Session = GUDUtils.get_engine_session(db_name)
feat = ClinVar.select_by_clinvarID(Session, None, "5000").first()
ClinVar.as_genomic_feature(feat)
