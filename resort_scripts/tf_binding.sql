CREATE TABLE `tf_binding_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `tf` varchar(25) NOT NULL,
  `score` float DEFAULT NULL,
  `peak` int(11) DEFAULT NULL,
  `regionID` int(10) unsigned NOT NULL,
  `sampleID` int(10) unsigned NOT NULL,
  `experimentID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`sampleID`,`experimentID`,`sourceID`,`tf`,`peak`),
  KEY `sampleID` (`sampleID`),
  KEY `experimentID` (`experimentID`),
  KEY `sourceID` (`sourceID`),
  KEY `ix_join` (`regionID`,`sampleID`,`experimentID`,`sourceID`),
  CONSTRAINT `tf_binding_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `tf_binding_sorted_ibfk_2` FOREIGN KEY (`sampleID`) REFERENCES `samples` (`uid`),
  CONSTRAINT `tf_binding_sorted_ibfk_3` FOREIGN KEY (`experimentID`) REFERENCES `experiments` (`uid`),
  CONSTRAINT `tf_binding_sorted_ibfk_4` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

INSERT INTO tf_binding_sorted
SELECT @rn:=@rn+1 AS uid, t1.tf, t1.score, t1.peak, t1.regionID, t1.sampleID, t1.experimentID, t1.sourceID
FROM (
  SELECT t.*
  from tf_binding t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;