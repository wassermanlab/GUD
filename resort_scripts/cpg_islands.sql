CREATE TABLE `cpg_islands_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `cpgs` int(11) NOT NULL,
  `gcs` int(11) NOT NULL,
  `percent_cpg` float NOT NULL,
  `percent_gc` float NOT NULL,
  `obsexp_ratio` float NOT NULL,
  `regionID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`sourceID`),
  KEY `sourceID` (`sourceID`),
  KEY `ix_join` (`regionID`,`sourceID`),
  CONSTRAINT `cpg_islands_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `cpg_islands_sorted_ibfk_2` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

INSERT INTO cpg_islands_sorted
SELECT @rn:=@rn+1 AS uid, t1.cpgs,t1.gcs, t1.percent_cpg, t1.percent_gc, t1.obsexp_ratio, t1.regionID, t1.sourceID 
FROM (
  SELECT t.*
  from cpg_islands t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;