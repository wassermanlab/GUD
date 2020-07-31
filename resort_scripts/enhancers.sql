CREATE TABLE `enhancers_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regionID` int(10) unsigned NOT NULL,
  `sampleID` int(10) unsigned NOT NULL,
  `experimentID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`sampleID`,`experimentID`,`sourceID`),
  KEY `sampleID` (`sampleID`),
  KEY `experimentID` (`experimentID`),
  KEY `sourceID` (`sourceID`),
  KEY `ix_join` (`regionID`,`sampleID`,`experimentID`,`sourceID`),
  CONSTRAINT `enhancers_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `enhancers_sorted_ibfk_2` FOREIGN KEY (`sampleID`) REFERENCES `samples` (`uid`),
  CONSTRAINT `enhancers_sorted_ibfk_3` FOREIGN KEY (`experimentID`) REFERENCES `experiments` (`uid`),
  CONSTRAINT `enhancers_sorted_ibfk_4` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

INSERT INTO enhancers_sorted
SELECT @rn:=@rn+1 AS uid, t1.regionID, t1.sampleID, t1.experimentID, t1.sourceID
FROM (
  SELECT t.*
  from enhancers t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;