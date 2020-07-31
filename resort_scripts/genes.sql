CREATE TABLE `genes_sorted` (
  `uid` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(75) NOT NULL,
  `gene_symbol` varchar(75) NOT NULL,
  `coding_start` int(10) unsigned NOT NULL,
  `coding_end` int(10) unsigned NOT NULL,
  `exon_starts` longblob NOT NULL,
  `exon_ends` longblob NOT NULL,
  `strand` char(1) DEFAULT NULL,
  `regionID` int(10) unsigned NOT NULL,
  `sourceID` int(10) unsigned NOT NULL,
  PRIMARY KEY (`uid`),
  UNIQUE KEY `regionID` (`regionID`,`name`,`sourceID`,`strand`),
  KEY `ix_join` (`sourceID`,`regionID`),
  KEY `ix_gene_symbol` (`gene_symbol`),
  KEY `ix_name` (`name`),
  CONSTRAINT `genes_sorted_ibfk_1` FOREIGN KEY (`regionID`) REFERENCES `regions` (`uid`),
  CONSTRAINT `genes_sorted_ibfk_2` FOREIGN KEY (`sourceID`) REFERENCES `sources` (`uid`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

INSERT INTO genes_sorted
SELECT @rn:=@rn+1 AS uid, t1.name, t1.gene_symbol, t1.coding_start, t1.coding_end, t1.exon_starts, t1.exon_ends, t1.strand, t1.regionID, t1.sourceID 
FROM (
  SELECT t.*
  from genes t, regions where t.regionID = regions.uid order by chrom, start, end limit 900000000
) t1, (SELECT @rn:=0) t2;