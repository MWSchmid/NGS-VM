# NGS-VM
Deploy a machine for NGS data analysis.

Note: I used this playbook to setup my NGS VMs on the ScienceCloud at the UZH. The instructions below are working on this specific infrastructure. Important: Ansible is not yet compatible with Python 3. Ubuntu 16.04 does not have Python 2.7 pre-installed anymore. If you want to use Ubuntu 16.04, you need to either write a [workaround](https://groups.google.com/forum/#!topic/ansible-project/DUKzTho3OCI) in Ansible, or modify the machine manually before running the playbook. Below I used Ubuntu 14.04 (server). The default-user is called "ubuntu" and there is no password (needs to be added on machines with a password).

Warning: I tested the VM-setup in June 2016 (I did the first version end of 2015). I had to update two packages (SPP and MULTOVL) and still need to check/test whether the updates have an effect on the functionality (phantompeakqualtools requires SPP). In general, some parts may become outdated in relatively short time... Another thing one could think of is to set the version number of the packages which are directly downloaded with pip.

While testing it again, building picard failed. Will be fixed in the next days.

## Requirements

1. Install dependencies:
```sh
sudo apt-get install git curl ansible ssh
```
2. Clone the repository:
```sh
git clone https://github.com/MWSchmid/NGS-VM
```
3. GitHub is not meant to store large files just for convenience. You therefore have to download all the packages yourself. This is only required once but it is also prone to "update-errors". The script I provided contains links to all the packages which can be quickly outdated (the errors will normally result in "tar: Unrecognized archive format"). You eventually need to change some of the links in future. However, this may also change the requirements in general.
```sh
cd NGS-VM
# download all tools
sh downloadAllVMpackages.sh
# move the packages to the right place
for ARCHIVE in bamtools bedtools bismark bowtie1 bowtie2 bwa fastqc fqtrim HTSeq multovl Rcount RSEM samtools soapAligner soapBuilder star subread trimGalore UCSC; do
mv repacked/$ARCHIVE.tar.gz roles/commonTools/files/
done
for ARCHIVE in CEAS GEM htsjdk jChIP MACS MAnorm PAPST PeakAnalyzer phantompeakqualtools picard SICER SPP trimmomatic; do
mv repacked/$ARCHIVE.tar.gz roles/additionalTools/files/
done
```

## Deploy a VM

* Launch a VM
* Replace the IP in [hosts](hosts) with the one from the VM. To deploy more VMs at once, just add the other IPs in this file as well.
* If not done yet, do:
```sh
ssh-add /path/to/my/private/key.pem
```
* Edit [setupVM.yml](setupVM.yml) if necessary (you can for example disable some of the roles) and start the deployment process:
```sh
ansible-playbook setupVM.yml -i hosts -e "uservar=ubuntu"
```
* To enable other users on the VM, log into the VM and remove some user specific files:
```sh
ssh ubuntu@172.23.xxx.xxx
sudo rm /var/lib/cloud/data/*
sudo rm -rf /var/lib/cloud/instances/*
exit
```

## Notes

All the files generated during the setup are located under "/home/setup". Most important binaries were moved to "/usr/local/bin" - except for the UCSC tools which can be found under "/home/setup/UCSC".

Finally, some [generic examples](genericExamples.md) for some tasks.



