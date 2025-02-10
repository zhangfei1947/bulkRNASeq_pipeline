#!/usr/bin/env python3
import os
import yaml
import sys
from pathlib import Path

class ConfigValidator:
    def __init__(self, config_path):
        self.config_path = config_path
        self.errors = []
        self.warnings = []
        self.config = None
        
    def load_config(self):
        try:
            with open(self.config_path) as f:
                self.config = yaml.safe_load(f)
        except Exception as e:
            self.errors.append(f"Failed to load config file: {str(e)}")
            return False
        return True
    
    def validate_structure(self):
        required_sections = ['samples', 'genome', 'diff_comparisons', 'venn_genes', 'enrichment']
        for section in required_sections:
            if section not in self.config:
                self.errors.append(f"Missing required section: [{section}]")
                
    def validate_samples(self):
        samples = self.config.get('samples', {})
        if not samples:
            self.errors.append("No samples defined in [samples] section")
            return
            
        valid_groups = set()
        for sample_name, sample_info in samples.items():
            # Check required fields
            required_fields = ['raw_dir', 'raw_base', 'group']
            for field in required_fields:
                if field not in sample_info:
                    self.errors.append(f"Sample [{sample_name}] missing required field: {field}")
                    
            # Check raw data directory
            if 'raw_dir' in sample_info:
                if not Path(sample_info['raw_dir']).exists():
                    self.warnings.append(f"Raw data directory not found for sample [{sample_name}]: {sample_info['raw_dir']}")
                    
            valid_groups.add(sample_info.get('group', ''))
            
        self.valid_groups = valid_groups - {''}
        
    def validate_genome(self):
        genome = self.config.get('genome', {})
        required_genome = ['index', 'annotation']
        for field in required_genome:
            if field not in genome:
                self.errors.append(f"Missing genome.{field} configuration")
                
        # Check index files
        if 'index' in genome:
            index_base = genome['index']
            for ext in ['.1.ht2', '.2.ht2']:
                if not Path(f"{index_base}{ext}").exists():
                    self.errors.append(f"Hisat2 index file not found: {index_base}{ext}")
                    
        # Check annotation file
        if 'annotation' in genome:
            anno_path = Path(genome['annotation'])
            if not anno_path.exists():
                self.errors.append(f"Annotation file not found: {anno_path}")
            elif anno_path.suffix not in ['.gtf', '.gff']:
                self.warnings.append(f"Unexpected annotation file format: {anno_path.suffix}")
                
    def validate_comparisons(self):
        comparisons = self.config.get('diff_comparisons', {})
        if not comparisons:
            self.errors.append("No differential comparisons defined")
            return
            
        for comp_name, comp_info in comparisons.items():
            required_fields = ['numerator', 'denominator']
            for field in required_fields:
                if field not in comp_info:
                    self.errors.append(f"Comparison [{comp_name}] missing field: {field}")
                    
            # Check group existence
            for field in ['numerator', 'denominator']:
                group = comp_info.get(field, '')
                if group and group not in self.valid_groups:
                    self.errors.append(f"Invalid group '{group}' in comparison [{comp_name}]. Valid groups: {self.valid_groups}")
                    
    def validate_venn(self):
        venn_genes = self.config.get('venn_genes', [])
        valid_comparisons = set(self.config['diff_comparisons'].keys())
        
        for venn_set in venn_genes:
            if not isinstance(venn_set, dict):
                self.errors.append("Invalid venn_genes format, should be list of dictionaries")
                continue
                
            for name, sources in venn_set.items():
                for source in sources:
                    if isinstance(source, str) and '_vs_' in source:
                        if source not in valid_comparisons:
                            self.errors.append(f"Venn reference to undefined comparison: {source}")
                    elif isinstance(source, dict):
                        pass  # Can add custom validation for filter criteria
    
    def validate_enrichment(self):
        enrichment = self.config.get('enrichment', {})
        valid_go_dbs = ['BP', 'MF', 'CC']
        valid_kegg_org = ['hsa', 'mmu', 'rno']  # Add more as needed
        
        # Validate GO
        go_conf = enrichment.get('GO', {})
        if go_conf.get('databases', []):
            for db in go_conf['databases']:
                if db not in valid_go_dbs:
                    self.errors.append(f"Invalid GO database: {db}. Valid options: {valid_go_dbs}")
                    
        # Validate KEGG
        kegg_conf = enrichment.get('KEGG', {})
        org_db = kegg_conf.get('org_db', '')
        if org_db and org_db not in valid_kegg_org:
            self.errors.append(f"Unsupported KEGG organism: {org_db}. Valid options: {valid_kegg_org}")
            
    def run_validation(self):
        if not self.load_config():
            return False
            
        self.validate_structure()
        self.validate_samples()
        self.validate_genome()
        self.validate_comparisons()
        self.validate_venn()
        self.validate_enrichment()
        
        return len(self.errors) == 0
    
    def print_report(self):
        print(f"\nValidation Report for {self.config_path}")
        print("="*50)
        
        if self.errors:
            print("\n[ERRORS] Need immediate attention:")
            for i, err in enumerate(self.errors, 1):
                print(f"{i}. {err}")
                
        if self.warnings:
            print("\n[WARNINGS] Should be checked:")
            for i, warn in enumerate(self.warnings, 1):
                print(f"{i}. {warn}")
                
        if not self.errors and not self.warnings:
            print("\nConfig file is valid! No issues found.")
            
        print("="*50)
        
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: config_validate.py <config.yaml>")
        sys.exit(1)
        
    validator = ConfigValidator(sys.argv[1])
    is_valid = validator.run_validation()
    validator.print_report()
    
    if validator.errors:
        sys.exit(1)
    else:
        print("Validation passed successfully")
        sys.exit(0)

