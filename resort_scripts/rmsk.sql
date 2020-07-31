CREATE TABLE `rmsk_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `score` int(11) NOT NULL,
  `name` varchar(75) NOT NULL,
  `class` varchar(75) NOT NULL,
  `family` varchar(75) NOT NULL,
  `strand` char(1) DEFAULT NULL,
  `regionID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`sourceID`,`name`,`strand`),
  KEY `sourceID` (`sourceID`),
  KEY `ix_family` (`family`),
  KEY `ix_class` (`class`),
  KEY `ix_join` (`regionID`,`sourceID`),
  CONSTRAINT `rmsk_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `rmsk_sorted_ibfk_2` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB AUTO_INCREMENT=5317292 DEFAULT CHARSET=utf8;

INSERT INTO rmsk_sorted
SELECT @rn:=@rn+1 AS uid, t1.score, t1.name, t1.class, t1.family, t1.strand, t1.regionID, t1.sourceID
FROM (
  SELECT t.*
  from rmsk t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;