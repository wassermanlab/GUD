CREATE TABLE `conservation_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `score` float DEFAULT NULL,
  `regionID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`sourceID`),
  KEY `ix_join` (`sourceID`,`regionID`),
  CONSTRAINT `conservation_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `conservation_sorted_ibfk_2` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

INSERT INTO conservation_sorted
SELECT @rn:=@rn+1 AS uid, t1.score, t1.regionID, t1.sourceID 
FROM (
  SELECT t.*
  from conservation t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;